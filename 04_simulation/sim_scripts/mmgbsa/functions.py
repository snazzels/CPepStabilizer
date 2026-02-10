import os
import sys
import re
import mdtraj as md
import numpy as np
import pandas as pd
from mdtraj.geometry import _geometry
from mdtraj.utils import ensure_type

gbsa_pattern   = r'DELTA TOTAL\s+(?P<Average>-?\d+.\d+)\s+(?P<Std>\d+.\d+)\s+(?P<SEM>\d+.\d+)'
decomp_pattern = '\\S{3,}\\s+(?P<idx>\\d+),[RL]\\s\\S{3}\\s+\\d+,(-?\\d+.\\d+e?-?[0-9]*,){15}(?P<energy>-?\\d+.\\d+e?-?[0-9]*)'

def read_output(path,pattern):
    with open(path,'r') as file:
        content = file.read()
        matches = re.findall(pattern,content,re.MULTILINE)
        if isinstance(matches[0], tuple):
            return [float(match) for match in matches[0]]
        else:
            return [float(match) for match in matches]

# MD simulation analysis
class ABL:
    def __init__(self,sim_path,name,runs=5):
        ## Initialize
        self.runs = runs
        self.sim_path = sim_path
        self.name = name
        self.prm = ''.join([sim_path,name,'/',name,'.prmtop'])

    def analysis(self,A_resid,B_resid,A_intres,B_intres):
        ## Setup arguments
        Ares_arg = ''.join([':',str(min(A_resid)),'-',str(max(A_resid))])
        Aintres_arg = ':'+','.join([str(resid) for resid in A_intres])
        Bintres_arg = ':'+','.join([str(resid) for resid in B_intres])
        Bres_arg = ''.join([':',str(min(B_resid)),'-',str(max(B_resid))])
        Lres_arg = ''.join([':',str(max(B_resid)+1)])
    
    ## Read MMGBSA output, including total ddG and energy decompositions
    def read_mmgbsa(self,gbsa_dir='../Step6_MMGBSA/'):
        
        # initialize
        self.AL_gbsa  = np.zeros(self.runs)
        self.BL_gbsa  = np.zeros(self.runs)
        
        self.AL_gbsa_decomp  = [np.array] * self.runs
        self.BL_gbsa_decomp  = [np.array] * self.runs
        
        # read value
        for run in range(1,self.runs+1):
            # A
            AL_gbsa_path   = ''.join([gbsa_dir,self.name,'/AL_output_',str(run),'.dat'])
            AL_decomp_path = ''.join([gbsa_dir,self.name,'/AL_decomp_',str(run),'.csv'])
            self.AL_gbsa[run-1]        = float(read_output(AL_gbsa_path,gbsa_pattern)[0][0])
            self.AL_gbsa_decomp[run-1] = pd.read_csv(AL_decomp_path,skiprows=8,header=None)[17].values[:-1]
            
            # B
            BL_gbsa_path = ''.join([gbsa_dir,self.name,'/BL_output_',str(run),'.dat'])
            BL_decomp_path = ''.join([gbsa_dir,self.name,'/BL_decomp_',str(run),'.csv'])
            self.BL_gbsa[run-1] = float(read_output(BL_gbsa_path,gbsa_pattern)[0][0])
            self.BL_gbsa_decomp[run-1] = pd.read_csv(BL_decomp_path,skiprows=8,header=None)[17].values[:-1]
            
        
      

   
