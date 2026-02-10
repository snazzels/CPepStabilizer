import jax
import jax.numpy as jnp
import numpy as np
from jax import vmap
from colabdesign.af.alphafold.common import residue_constants

from functions_design.constants import VAN_DER_WAALS_RADII

def generate_sphere_points(n_sphere_points=960):
    """Generate evenly distributed points on a unit sphere using the Golden Section Spiral method."""
    indices = jnp.arange(0, n_sphere_points, dtype=jnp.float32)
    phi = jnp.pi * (3.0 - jnp.sqrt(5.0))  # Golden angle
    y = 1 - (indices / float(n_sphere_points - 1)) * 2  # y from 1 to -1
    radius = jnp.sqrt(1 - y**2)  # Corresponding xz radius
    theta = phi * indices  # Spiral angle

    x = radius * jnp.cos(theta)
    z = radius * jnp.sin(theta)
    
    return jnp.stack([x, y, z], axis=-1)  # (n_sphere_points, 3)

def calculate_sasa(atom_positions, van_der_waals_radii=VAN_DER_WAALS_RADII, probe_radius=1.4, n_sphere_points=960):
    """Compute Solvent Accessible Surface Area (SASA) using the Shrake-Rupley method."""
    
    sphere_points = generate_sphere_points(n_sphere_points)
    
    # Flatten atom positions for distance checking
    all_atoms = jnp.concatenate([p for p in atom_positions.values()], axis=0)
    
    # Extract radii for corresponding atoms
    filtered_keys = [key for key in atom_positions if key in van_der_waals_radii]
    all_radii = jnp.array([van_der_waals_radii[atom] for atom in filtered_keys for _ in atom_positions[atom]])

    def compute_accessible_surface(atom_pos, atom_radius):
        """Compute SASA for a single atom, skipping atoms with zero coordinates."""
        is_zero = jnp.all(atom_pos == 0)  # JAX-friendly condition

        # Generate sphere points around this atom
        expanded_sphere_points = sphere_points * (atom_radius + probe_radius)
        test_positions = atom_pos + expanded_sphere_points  # Shift sphere points

        # Compute distances to all other atoms
        distances = jnp.linalg.norm(test_positions[:, None, :] - all_atoms[None, :, :], axis=-1)

        # Check occlusion (True if distance < other atom's vdW radius) + buffer to not self exclude
        is_occluded = jnp.any((distances + 1e-5) < (all_radii[None, :] + probe_radius), axis=1)

        # Accessible points (not occluded)
        accessible_points = jnp.sum(~is_occluded)
        
        # Compute SASA contribution
        surface_area_per_point = 4 * jnp.pi * (atom_radius + probe_radius) ** 2 / n_sphere_points

        sasa_per_atom = jnp.where(is_zero, 0.0, accessible_points * surface_area_per_point)

        return sasa_per_atom

    # Vectorized computation over all atoms
    sasa_values = vmap(compute_accessible_surface)(all_atoms, all_radii)
    
    return jnp.sum(sasa_values)  # Total SASA

@jax.jit
def sasa_loss(inputs, outputs):
    """
    Compute normalized SASA loss for peptide binding to a protein-protein interface.
    
    The loss function has three components:
    1. BSA difference term - minimizes difference in binding between proteins A and B
    2. BSA magnitude term - guides BSA toward optimal values
    3. Minimum BSA threshold - ensures sufficient binding
    
    Returns:
    - Dictionary containing loss values
    """
    # Default parameters
    w1 = 1.0  # Weight for BSA difference
    w2 = 0.5  # Weight for total BSA
    w3 = 2.0  # Weight for minimum BSA threshold
    min_bsa_threshold = 500.0  # Minimum BSA in Å²
    
    final_positions = outputs["structure_module"]["final_atom_positions"]
    
    # Convert final positions into a dictionary format based on residue_constants
    chain_A_atoms = {atom: final_positions[:, idx][:276] for atom, idx in residue_constants.atom_order.items() if idx < final_positions.shape[1]} 
    chain_B_atoms = {atom: final_positions[:, idx][276:-8] for atom, idx in residue_constants.atom_order.items() if idx < final_positions.shape[1]} 
    peptide_atoms = {atom: final_positions[:, idx][-8:] for atom, idx in residue_constants.atom_order.items() if idx < final_positions.shape[1]}
    
    # Combine atoms for interface calculations
    chain_A_peptide_atoms = {}
    for atom in chain_A_atoms:
        chain_A_peptide_atoms[atom] = jnp.concatenate([chain_A_atoms[atom], peptide_atoms[atom]], axis=0) 
    chain_B_peptide_atoms = {}
    for atom in chain_B_atoms:
        chain_B_peptide_atoms[atom] = jnp.concatenate([chain_B_atoms[atom], peptide_atoms[atom]], axis=0) 
  
    # Calculate SASA for individual components and complexes
    sasa_chain_A = calculate_sasa(chain_A_atoms)
    sasa_chain_B = calculate_sasa(chain_B_atoms)
    sasa_peptide = calculate_sasa(peptide_atoms)
    sasa_chain_A_peptide = calculate_sasa(chain_A_peptide_atoms)
    sasa_chain_B_peptide = calculate_sasa(chain_B_peptide_atoms)
    
    # Calculate buried surface area (BSA) for each protein-peptide interface
    bsa_a = sasa_chain_A + sasa_peptide - sasa_chain_A_peptide
    bsa_b = sasa_chain_B + sasa_peptide - sasa_chain_B_peptide
    
    # 1. BSA difference term - penalize unequal binding
    max_possible_diff = sasa_peptide  # Maximum possible difference can't exceed total peptide SASA
    bsa_difference = jnp.abs(bsa_a - bsa_b) / max_possible_diff
    
    # 2. Total BSA term - we want to maximize total binding (negative as we want to maximize)
    total_bsa = bsa_a + bsa_b
    max_possible_bsa = 2.0 * sasa_peptide  # Theoretical maximum if peptide was completely buried
    bsa_magnitude_term = 1.0 - (total_bsa / max_possible_bsa)  # Normalized to [0,1]
    
    # 3. Minimum BSA penalty - ensure each interface has sufficient binding
    min_bsa = jnp.minimum(bsa_a, bsa_b)
    min_bsa_penalty = jnp.maximum(0.0, 1.0 - (min_bsa / min_bsa_threshold))
    
    # Combine the terms with weights
    normalized_loss = (w1 * bsa_difference + 
                      w2 * bsa_magnitude_term + 
                      w3 * min_bsa_penalty) / (w1 + w2 + w3)
    
    # Return a dictionary with detailed loss components for analysis
    return {"bsa_loss": normalized_loss}

