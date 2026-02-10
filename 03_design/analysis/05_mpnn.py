#!/usr/bin/env python3


import os
import re
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import softmax
from scipy.stats import pearsonr, spearmanr

# Import necessary libraries from colabdesign
import colabdesign
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.shared.protein import pdb_to_string
from colabdesign.af.alphafold.common import residue_constants

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Amino acid order and reverse mapping
aa_order = residue_constants.restype_order
order_aa = {b:a for a,b in aa_order.items()}

# Similar residues dictionary for adjusted sequence identity
SIMILAR_RESIDUES = {
    'F': ['W', 'Y'], 'W': ['F', 'Y'], 'Y': ['F', 'W'],  
    'L': ['I', 'V'], 'I': ['L', 'V'], 'V': ['L', 'I'],
    'A': ['G'], 'G': ['A'], 'S': ['T'], 'T': ['S'],
    'N': ['Q'], 'Q': ['N'], 'K': ['R'], 'R': ['K'],
    'D': ['E'], 'E': ['D']
}

def analyze_pssm_sequence_correlation(pssm, sequences, n_residues=10, output_file=None):
    """
    Analyze correlation between PSSM probabilities and actual sampled sequences.
    
    Parameters:
    - pssm: numpy array of shape (positions, amino_acids), softmaxed probabilities
    - sequences: list of sequence strings (the 32 sampled sequences)
    - n_residues: number of residues from the end to analyze
    - output_file: optional file to write detailed analysis
    
    Returns:
    - correlation_stats: dictionary with correlation statistics
    """
    # Amino acid mapping
    #amino_acids = list("ACDEFGHIKLMNPQRSTVWYX")
    amino_acids = list("ARNDCQEGHILKMFPSTWYVX")

    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
    
    # Get last n_residues positions
    pssm_subset = pssm[-n_residues:]
    
    # Calculate frequency of each amino acid at each position from sequences
    sequence_frequencies = np.zeros((n_residues, len(amino_acids)))
    
    for seq in sequences:
        seq_subset = seq[-n_residues:]  # Last n_residues of this sequence
        for pos_idx, aa in enumerate(seq_subset):
            if aa in aa_to_idx:
                sequence_frequencies[pos_idx, aa_to_idx[aa]] += 1
    
    # Normalize to get frequencies (probabilities)
    sequence_frequencies = sequence_frequencies / len(sequences)
    
    # Calculate correlations
    correlations = []
    position_stats = []
    
    for pos_idx in range(n_residues):
        pssm_probs = pssm_subset[pos_idx]
        seq_freqs = sequence_frequencies[pos_idx]
        
        # Calculate Pearson correlation
        pearson_corr, pearson_p = pearsonr(pssm_probs, seq_freqs)
        spearman_corr, spearman_p = spearmanr(pssm_probs, seq_freqs)
        
        correlations.append({
            'position': len(pssm) - n_residues + pos_idx,
            'pearson_corr': pearson_corr,
            'pearson_p': pearson_p,
            'spearman_corr': spearman_corr,
            'spearman_p': spearman_p
        })
        
        # Find top amino acids in PSSM vs sequences
        pssm_top_indices = np.argsort(pssm_probs)[-5:][::-1]  # Top 5
        seq_top_indices = np.argsort(seq_freqs)[-5:][::-1]    # Top 5
        
        pssm_top_aa = [(amino_acids[i], pssm_probs[i]) for i in pssm_top_indices]
        seq_top_aa = [(amino_acids[i], seq_freqs[i]) for i in seq_top_indices]
        
        position_stats.append({
            'position': len(pssm) - n_residues + pos_idx,
            'pssm_top': pssm_top_aa,
            'sequence_top': seq_top_aa
        })
    
    # Overall statistics
    avg_pearson = np.mean([c['pearson_corr'] for c in correlations])
    avg_spearman = np.mean([c['spearman_corr'] for c in correlations])
    
    correlation_stats = {
        'avg_pearson_correlation': avg_pearson,
        'avg_spearman_correlation': avg_spearman,
        'position_correlations': correlations,
        'position_details': position_stats,
        'overall_pssm_entropy': np.mean(-np.sum(pssm_subset * np.log(pssm_subset + 1e-9), axis=1)),
        'overall_sequence_entropy': np.mean(-np.sum(sequence_frequencies * np.log(sequence_frequencies + 1e-9), axis=1))
    }
    
    # Write detailed analysis to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write("PSSM vs Sequence Correlation Analysis\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Average Pearson Correlation: {avg_pearson:.4f}\n")
            f.write(f"Average Spearman Correlation: {avg_spearman:.4f}\n")
            f.write(f"PSSM Entropy: {correlation_stats['overall_pssm_entropy']:.4f}\n")
            f.write(f"Sequence Entropy: {correlation_stats['overall_sequence_entropy']:.4f}\n\n")
            
            f.write("Position-by-Position Analysis:\n")
            f.write("-" * 30 + "\n")
            
            for i, (corr, detail) in enumerate(zip(correlations, position_stats)):
                f.write(f"\nPosition {detail['position']}:\n")
                f.write(f"  Pearson r: {corr['pearson_corr']:.4f} (p={corr['pearson_p']:.4f})\n")
                f.write(f"  Spearman r: {corr['spearman_corr']:.4f} (p={corr['spearman_p']:.4f})\n")
                
                f.write("  Top 5 PSSM probabilities: ")
                f.write(", ".join([f"{aa}:{prob:.3f}" for aa, prob in detail['pssm_top']]))
                f.write("\n")
                
                f.write("  Top 5 sequence frequencies: ")
                f.write(", ".join([f"{aa}:{freq:.3f}" for aa, freq in detail['sequence_top']]))
                f.write("\n")
                
                # Check overlap in top predictions
                pssm_top_aa_only = [aa for aa, _ in detail['pssm_top']]
                seq_top_aa_only = [aa for aa, _ in detail['sequence_top']]
                overlap = set(pssm_top_aa_only) & set(seq_top_aa_only)
                f.write(f"  Top-5 overlap: {len(overlap)}/5 amino acids ({', '.join(overlap)})\n")
    
    return correlation_stats

def create_correlation_plot(pssm, sequences, n_residues=10, output_file=None):
    """
    Create a scatter plot showing PSSM probabilities vs sequence frequencies.
    """
    #amino_acids = list("ACDEFGHIKLMNPQRSTVWYX")
    amino_acids = list("ARNDCQEGHILKMFPSTWYVX")

    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
    
    # Get data
    pssm_subset = pssm[-n_residues:]
    
    # Calculate sequence frequencies
    sequence_frequencies = np.zeros((n_residues, len(amino_acids)))
    for seq in sequences:
        seq_subset = seq[-n_residues:]
        for pos_idx, aa in enumerate(seq_subset):
            if aa in aa_to_idx:
                sequence_frequencies[pos_idx, aa_to_idx[aa]] += 1
    sequence_frequencies = sequence_frequencies / len(sequences)
    
    # Create scatter plot
    plt.figure(figsize=(10, 8))
    
    # Flatten arrays for plotting
    pssm_flat = pssm_subset.flatten()
    freq_flat = sequence_frequencies.flatten()
    
    # Color by amino acid type
    colors = []
    labels = []
    for pos in range(n_residues):
        for aa_idx, aa in enumerate(amino_acids):
            colors.append(pos)
            labels.append(f"Pos{len(pssm)-n_residues+pos}_{aa}")
    
    scatter = plt.scatter(pssm_flat, freq_flat, c=colors, alpha=0.6, cmap='viridis')
    plt.colorbar(scatter, label='Position')
    
    # Add diagonal line for perfect correlation
    max_val = max(np.max(pssm_flat), np.max(freq_flat))
    plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='Perfect correlation')
    
    plt.xlabel('PSSM Probability')
    plt.ylabel('Sequence Frequency')
    plt.title(f'PSSM Probabilities vs Sequence Frequencies (Last {n_residues} residues)')
    plt.legend()
    
    # Calculate and display overall correlation
    overall_pearson, _ = pearsonr(pssm_flat, freq_flat)
    plt.text(0.05, 0.95, f'Overall Pearson r = {overall_pearson:.3f}', 
             transform=plt.gca().transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return plt.gcf()

def extract_top_probabilities(pssm, n_residues=0, top_k=5):
    """
    Extract the top k probabilities for the last n residues
    
    Parameters:
    - pssm: numpy array of shape (positions, amino_acids)
    - n_residues: number of residues from the end to analyze
    - top_k: number of top probabilities to keep for each position
    
    Returns:
    - top_probs: dictionary with position keys and list of (amino_acid, probability) tuples
    """
    # Get the last n_residues positions
    last_positions = pssm[-n_residues:]
    
    # Dictionary to store results
    top_probs = {}
    
    # Amino acid list (assuming standard 20 amino acids + X)
    #amino_acids = list("ACDEFGHIKLMNPQRSTVWYX")
    amino_acids = list("ARNDCQEGHILKMFPSTWYVX")

    
    # Process each position
    for i, pos_data in enumerate(last_positions):
        position_idx = len(pssm) - n_residues + i
        
        # Create pairs of (amino acid, probability)
        aa_prob_pairs = [(aa, prob) for aa, prob in zip(amino_acids, pos_data)]
        
        # Sort by probability in descending order and keep top k
        top_k_pairs = sorted(aa_prob_pairs, key=lambda x: x[1], reverse=True)[:top_k]
        
        # Store in dictionary
        top_probs[f"Pos_{position_idx}"] = top_k_pairs
    
    return top_probs

def write_top_probabilities_to_file(top_probs, output_file):
    """
    Write the top probabilities to a file in a formatted table.
    
    Parameters:
    - top_probs: Dictionary with position keys and list of (amino_acid, probability) tuples
    - output_file: Path to the output file
    """
    with open(output_file, "w") as f:
        headers = ["Position"] + [f"Rank {i+1}" for i in range(len(list(top_probs.values())[0]))]
        f.write(" | ".join(headers) + "\n")
        f.write("-" * (len(headers) * 10) + "\n")

        for position, probs in top_probs.items():
            row = [position] + [f"{aa} ({prob:.3f})" for aa, prob in probs]
            f.write(" | ".join(row) + "\n")

def calculate_sequence_identity(sequences, reference_seq, similar_residues=SIMILAR_RESIDUES):
    """
    Calculate exact and adjusted sequence identity.
    
    Parameters:
    - sequences: list or series of sequences to compare
    - reference_seq: reference sequence to compare against
    - similar_residues: dictionary of similar amino acids
    
    Returns:
    - tuple of (exact_identity, adjusted_identity)
    """
    n_res = len(reference_seq)
    
    # Calculate exact and adjusted matches
    exact_matches = sequences.apply(lambda seq: sum(1 for a, b in zip(seq, reference_seq) if a == b))
    adjusted_matches = sequences.apply(lambda seq: sum(1 for a, b in zip(seq, reference_seq) 
                                                       if a == b or 
                                                       (b in similar_residues and a in similar_residues[b])))
    
    exact_identity = exact_matches / n_res
    adjusted_identity = adjusted_matches / n_res
    
    return exact_identity, adjusted_identity


def calculate_perplexity(pssm, reference_seq, n_residues=10):
    """
    Calculate perplexity of reference sequence vs PSSM distribution
    for the last n residues.

    Parameters:
    - pssm: numpy array of shape (positions, amino_acids), already softmaxed
    - reference_seq: full reference sequence (string)
    - n_residues: number of residues from the end to evaluate

    Returns:
    - perplexity (float)
    """
    # Amino acid order used in pssm
    #amino_acids = list("ACDEFGHIKLMNPQRSTVWYX")
    amino_acids = list("ARNDCQEGHILKMFPSTWYVX")

    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}

    # Select last n residues and corresponding pssm slice
    ref_subseq = reference_seq[-n_residues:]
    pssm_sub = pssm[-n_residues:]

    # Collect probabilities assigned to the correct residues
    probs = []
    for aa, pos_probs in zip(ref_subseq, pssm_sub):
        if aa not in aa_to_idx:
            continue  # skip nonstandard residues
        probs.append(pos_probs[aa_to_idx[aa]])

    probs = np.array(probs)

    # Cross-entropy and perplexity
    log_likelihood = np.mean(np.log(probs + 1e-9))  # avoid log(0)
    perplexity = np.exp(-log_likelihood)

    return perplexity


def calculate_adjusted_perplexity(pssm, reference_seq, n_residues=10, similar_residues=SIMILAR_RESIDUES):
    """
    Calculate adjusted perplexity for the last n residues, 
    where probabilities of similar residues are also counted.
    """
    #amino_acids = list("ACDEFGHIKLMNPQRSTVWYX")
    amino_acids = list("ARNDCQEGHILKMFPSTWYVX")
    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}

    ref_subseq = reference_seq[-n_residues:]
    pssm_sub = pssm[-n_residues:]

    probs = []
    for aa, pos_probs in zip(ref_subseq, pssm_sub):
        if aa not in aa_to_idx:
            continue

        # Start with probability of exact match
        p = pos_probs[aa_to_idx[aa]]

        # Add probabilities of similar residues
        if aa in similar_residues:
            for sim_aa in similar_residues[aa]:
                if sim_aa in aa_to_idx:
                    p += pos_probs[aa_to_idx[sim_aa]]

        probs.append(p)

    probs = np.array(probs)
    log_likelihood = np.mean(np.log(probs + 1e-9))  # avoid log(0)
    adj_perplexity = np.exp(-log_likelihood)

    return adj_perplexity


def main(
    pdb_dir="cleaned_top/",
    output_dir="output",
    model_name="v_48_020",
    chains="A,B,C",
    homooligomer=False,
    fix_pos="A,B",
    inverse=False,
    rm_aa=None,
    num_seqs=32,
    sampling_temp=0.1,
    n_res_identity=14
):
    """
    Main function to run ProteinMPNN design and sequence identity analysis
    
    Parameters:
    - Various parameters for ProteinMPNN design and sequence identity calculation
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Prepare ProteinMPNN model
    mpnn_model = mk_mpnn_model(model_name)

    # Prepare results storage
    identity_data = []
    correlation_data = []

    # Process each PDB file
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    for pdb_path in pdb_files:
        pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
        print(f"\nProcessing {pdb_name}...")
        
        # Prepare MPNN inputs
        mpnn_model.prep_inputs(
            pdb_filename=pdb_path,
            chain=chains, 
            homooligomer=homooligomer,
            fix_pos=fix_pos, 
            inverse=inverse,
            rm_aa=rm_aa, 
            verbose=True
        )
        
        # Identify non-fixed sequence
        L = sum(mpnn_model._lengths)
        non_fixed_mask = np.ones(L, dtype=bool)
        if "fix_pos" in mpnn_model._inputs:
            non_fixed_mask[mpnn_model._inputs["fix_pos"]] = False
        non_fixed_seq = "".join([order_aa.get(a, "H") for a, keep in zip(mpnn_model._inputs["S"], non_fixed_mask) if keep])
        
        # Sample sequences
        out = mpnn_model.sample(
            num=num_seqs//32, 
            batch=32,
            temperature=sampling_temp,
            rescore=homooligomer
        )
        
        # Prepare output filenames
        fasta_filename = os.path.join(output_dir, f"{pdb_name}_design.fasta")
        csv_filename = os.path.join(output_dir, f"{pdb_name}.csv")
        
        # Write FASTA file
        with open(fasta_filename, "w") as fasta:
            for n in range(num_seqs):
                line = f'>score:{out["score"][n]:.3f}_seqid:{out["seqid"][n]:.3f}\n{out["seq"][n]}'
                fasta.write(line + "\n")
        
        # Prepare DataFrame
        labels = ["score", "seqid", "seq", "non_fixed_seq"]
        data = [[out["score"][n], out["seqid"][n], out["seq"][n], non_fixed_seq] for n in range(num_seqs)]
        df = pd.DataFrame(data, columns=labels)
        df.to_csv(csv_filename, index=False)
        
        # Calculate PSSM and save
        ar_mask = np.zeros((L,L))
        logits = mpnn_model.score(ar_mask=ar_mask)["logits"]
        pssm = softmax(logits, -1)
        np.savetxt(os.path.join(output_dir, f"{pdb_name}_pssm.txt"), pssm)

        # ===== NEW: PSSM-SEQUENCE CORRELATION ANALYSIS =====
        print("Analyzing PSSM-sequence correlations...")
        
        # Analyze correlation between PSSM and actual sequences
        correlation_analysis_file = os.path.join(output_dir, f"{pdb_name}_correlation_analysis.txt")
        correlation_stats = analyze_pssm_sequence_correlation(
            pssm, 
            out["seq"], 
            n_residues=n_res_identity,
            output_file=correlation_analysis_file
        )
        
        # Create correlation plot
        correlation_plot_file = os.path.join(output_dir, f"{pdb_name}_correlation_plot.png")
        correlation_fig = create_correlation_plot(
            pssm, 
            out["seq"], 
            n_residues=n_res_identity,
            output_file=correlation_plot_file
        )
        plt.close(correlation_fig)  # Close to free memory
        
        # Store correlation results
        correlation_data.append([
            pdb_name,
            correlation_stats['avg_pearson_correlation'],
            correlation_stats['avg_spearman_correlation'],
            correlation_stats['overall_pssm_entropy'],
            correlation_stats['overall_sequence_entropy']
        ])
        
        print(f"  Average Pearson correlation: {correlation_stats['avg_pearson_correlation']:.4f}")
        print(f"  Average Spearman correlation: {correlation_stats['avg_spearman_correlation']:.4f}")
        # ===== END NEW SECTION =====

        # Compute perplexity for last n residues
        perplexity = calculate_perplexity(pssm, non_fixed_seq, n_residues=n_res_identity)
        adj_perplexity = calculate_adjusted_perplexity(pssm, non_fixed_seq, n_residues=n_res_identity)

        print(f"  Perplexity (last {n_res_identity} residues): {perplexity:.3f}")
        print(f"  Adjusted perplexity: {adj_perplexity:.3f}")
        
        # Extract and write top probabilities
        top_probs = extract_top_probabilities(pssm, n_residues=n_res_identity)
        write_top_probabilities_to_file(top_probs, os.path.join(output_dir, f"{pdb_name}_top.txt"))
        
        # Extract sequences for last n_res_identity residues
        sequences = df['seq'].astype(str).str[-n_res_identity:]
        non_fixed_seq_subset = non_fixed_seq[-n_res_identity:]
        
        # Calculate sequence identities
        exact_identity, adjusted_identity = calculate_sequence_identity(sequences, non_fixed_seq_subset)
        
        # Store identity results
        avg_exact_identity = exact_identity.mean()
        avg_adjusted_identity = adjusted_identity.mean()
        
        identity_data.append([
            pdb_name, 
            avg_exact_identity, 
            avg_adjusted_identity,
            perplexity,
            adj_perplexity
        ])
        
        print(f"  Results saved to {csv_filename} and {fasta_filename}")
    
    # Create identity results DataFrame
    identity_df = pd.DataFrame(
        identity_data, 
        columns=['File', 'Avg Exact Identity', 'Avg Adjusted Identity', 'Perplexity', 'Adj. Perplexity']
    )
    
    # Create correlation results DataFrame
    correlation_df = pd.DataFrame(
        correlation_data,
        columns=['File', 'Avg Pearson Correlation', 'Avg Spearman Correlation', 'PSSM Entropy', 'Sequence Entropy']
    )
    
    # Sort and save results
    identity_df_sorted = identity_df.sort_values(by='Avg Exact Identity', ascending=False)
    identity_csv_path = os.path.join(output_dir, "sequence_identity_results.csv")
    identity_df_sorted.to_csv(identity_csv_path, index=False)
    
    correlation_df_sorted = correlation_df.sort_values(by='Avg Pearson Correlation', ascending=False)
    correlation_csv_path = os.path.join(output_dir, "correlation_results.csv")
    correlation_df_sorted.to_csv(correlation_csv_path, index=False)
    
    print(f"\nSaved sequence identity results to: {identity_csv_path}")
    print(f"Saved correlation analysis results to: {correlation_csv_path}")
    
    # Print summary statistics
    print(f"\nSUMMARY:")
    print(f"Average Pearson correlation across all files: {correlation_df['Avg Pearson Correlation'].mean():.4f}")
    print(f"Average Spearman correlation across all files: {correlation_df['Avg Spearman Correlation'].mean():.4f}")


if __name__ == "__main__":
    # Example usage with default parameters
    main()
