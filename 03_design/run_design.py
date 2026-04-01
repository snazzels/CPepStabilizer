import os
import re
import glob
import mdtraj as md
import io
import argparse
import shutil
import sys
from pathlib import Path
import datetime
import yaml

from colabdesign import mk_afdesign_model, clear_mem
from colabdesign.shared.utils import copy_dict
from colabdesign.af.alphafold.common import protein, residue_constants
import jax
import jax.numpy as jnp

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import softmax

# Custom loss functions (cyclic offset, COM distance, cis/trans penalty, BSA loss).
# Currently only add_cyclic_offset is active; the others are available for future use.
from functions_design.loss_functions import *
from functions_design.sasa_functions import *

# Load configuration
_repo_root = Path(__file__).resolve().parent.parent
with open(_repo_root / "config.yaml") as _f:
    _config = yaml.safe_load(_f)

class LogRedirector:
    """Redirect stdout and stderr to a file while still printing to console."""
    
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log_file = open(log_file, 'w')
    
    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush()
    
    def flush(self):
        self.terminal.flush()
        self.log_file.flush()
    
    def close(self):
        self.log_file.close()

def run_design(run_index, output_dir):
    run_dir = Path(output_dir) / f"run_{run_index:03d}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    run_dir.mkdir(exist_ok=True, parents=True)
    
    # Set up logging
    log_file = run_dir / "run.log"
    original_stdout = sys.stdout
    log_redirector = LogRedirector(log_file)
    sys.stdout = log_redirector
    
    try:
        print(f"\n\n{'='*80}")
        print(f"Starting run {run_index} in directory: {run_dir}")
        print(f"{'='*80}\n")
        
        # Load paths and parameters from config
        params_dir = _config["paths"]["af2_params_dir"]
        pdb = str(_repo_root / _config["paths"]["target_pdb"])

        design_cfg = _config["design"]
        target_chain = design_cfg["target_chain"]
        target_hotspot = design_cfg["target_hotspot"]
        target_hotspot_A = design_cfg["target_hotspot_A"]
        target_hotspot_B = design_cfg["target_hotspot_B"]

        if target_hotspot == "": target_hotspot = None
        target_flexible = False

        # Set binder length and initial AA sequence. set binder_seq = "" for random init. seq.
        binder_len = design_cfg["binder_len"]
        binder_seq = ""
        if len(binder_seq) > 0:
            binder_len = len(binder_seq)
        else:
            binder_seq = None

        # Offset to get cyclic and not linear peptide
        cyclic_offset = design_cfg["cyclic_offset"]

        # Bugfixed version of cyclic offset (True)
        bugfix = design_cfg["bugfix"]

        # PeptideLoss provides cis_loss and com_loss; instantiated but not currently
        # added to the model (available for future use via model.add_losses(...)).
        pep_loss = PeptideLoss(n_pep_res=binder_len, hotspot_res=target_hotspot, bound=100)

        # AF2 parameters (MCMC opt. uses only 1 model)
        use_multimer = design_cfg["use_multimer"]
        num_recycles = design_cfg["num_recycles"]
        num_models = design_cfg["num_models"]

        # Model parameters
        x = {"pdb_filename":pdb,
             "chain":target_chain,
             "binder_len":binder_len,
             "hotspot":target_hotspot,
             "hotspot_1":target_hotspot_A,
             "hotspot_2":target_hotspot_B,
             "use_multimer":use_multimer,
             "rm_target_seq":target_flexible}

        # Set model parameters and losses
        clear_mem()
        model = mk_afdesign_model(
            protocol="binder",
            use_multimer=x["use_multimer"],
            num_recycles=num_recycles,
            recycle_mode="sample", data_dir=params_dir)
        model.prep_inputs(**x, ignore_missing=False)
        model.opt["weights"]["plddt"] = 1.0
        print("weights", model.opt["weights"])
        print("target length:", model._target_len)
        print("binder length:", model._binder_len)
        binder_len = model._binder_len

        # Add cyclic offset if set True
        if cyclic_offset:
            if bugfix:
                print("Set bug fixed cyclic peptide complex offset. The cyclic peptide binder will be hallucinated.")
                add_cyclic_offset(model, bug_fix=True)
            else:
                print("Set not bug fixed cyclic peptide complex offset. The cyclic peptide binder will be hallucinated.")
                add_cyclic_offset(model, bug_fix=False)
        else:
            print("Don't set cyclic offset. The linear peptide binder will be hallucionated.")

        # Set optimizer for sequence optimization
        optimizer = design_cfg["optimizer"]

        # Advanced GD settings
        GD_method = "sgd"
        learning_rate = design_cfg["learning_rate"]
        norm_seq_grad = True
        dropout = True

        model.restart(seq=binder_seq)
        model.set_optimizer(optimizer=GD_method,
                            learning_rate=learning_rate,
                            norm_seq_grad=norm_seq_grad)
        models = model._model_names[:num_models]

        flags = {"num_recycles":num_recycles,
                 "models":models,
                 "dropout":dropout}

        if optimizer == "3stage":
            model.design_3stage(120, 60, 10, **flags)
            pssm = softmax(model._tmp["seq_logits"],-1)

        if optimizer == "my_pssm_semigreedy":
            model.my_design_pssm_semigreedy(120, 128, **flags)
            pssm = softmax(model._tmp["seq_logits"],1)

        if optimizer == "my_semigreedy":
            model.my_design_pssm_semigreedy(0, 256, **flags)
            pssm = None

        if optimizer == "pssm_semigreedy":
            model.design_pssm_semigreedy(80, 120, ramp_recycles=False, ramp_models=False, **flags)
            pssm = softmax(model._tmp["seq_logits"],1)

        if optimizer == "semigreedy":
            model.design_pssm_semigreedy(0, 256, ramp_recycles=False, ramp_models=False, **flags)
            pssm = None

        if optimizer == "mcmc":
            model._design_mcmc(steps=1000, mutation_rate=1, T_init=0.01,half_life=200, **flags)

        if optimizer == "pssm":
            model.design_logits(120, e_soft=1.0, num_models=1, ramp_recycles=True, **flags)
            model.design_soft(32, num_models=1, **flags)
            flags.update({"dropout":False,"save_best":True})
            model.design_soft(10, num_models=num_models, **flags)
            pssm = softmax(model.aux["seq"]["logits"],-1)

        # Save output to run directory instead of current directory
        output_pdb = run_dir / f"{model.protocol}.pdb"
        model.save_pdb(str(output_pdb))
        
        print(f"\nCompleted run {run_index}. Output saved to: {output_pdb}\n")
        
        # Save trajectory files if they exist
        trajectory_dir = Path("trajectory")
        if trajectory_dir.exists():
            dest_trajectory = run_dir / "trajectory"
            dest_trajectory.mkdir(exist_ok=True)
            for file in trajectory_dir.glob("*"):
                shutil.copy(file, dest_trajectory / file.name)
            
        return output_pdb
        
    except Exception as e:
        print(f"Error in run {run_index}: {str(e)}")
        return None
    finally:
        # Restore original stdout and close log file
        sys.stdout = original_stdout
        log_redirector.close()

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run peptide design script multiple times')
    parser.add_argument('-n', '--num_runs', type=int, default=5, 
                        help='Number of runs to perform (default: 5)')
    parser.add_argument('-o', '--output_dir', type=str, default='design_runs',
                        help='Base output directory (default: design_runs)')
    args = parser.parse_args()
    
    # Create base output directory
    base_dir = Path(args.output_dir)
    base_dir.mkdir(exist_ok=True, parents=True)
    
    results = []
    
    # Run the design process n times
    for i in range(1, args.num_runs + 1):
        output_pdb = run_design(i, base_dir)
        if output_pdb:
            results.append(output_pdb)
    
    print(f"\nAll {args.num_runs} runs completed. Results are in {base_dir}")
    print("Generated PDB files:")
    for pdb in results:
        print(f"  - {pdb}")

if __name__ == "__main__":
    main()
