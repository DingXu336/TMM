#!/usr/bin/env python3
"""
2D Jones Dispersion Demo - Reproducing original Jones_disperion.py functionality

This script creates 2D dispersion plots showing how optical properties vary with 
energy (eV) and q-vector (μm⁻¹), including:
- Reflectivity (Rs, Rp)
- Stokes parameters (S0, S1, S2, S3)  
- Polarization angles (phi, chi)

Author: Ding Xu (Physical Chemistry Researcher)
Based on original Jones_disperion.py
"""

import numpy as np
import matplotlib.pyplot as plt
from tmm.materials import get_builtin_material, Material
from tmm.utils.constants import Constants

# Physical constants 
c = Constants.C    # m/s
mu0 = Constants.MU0   # H/m
eps0 = Constants.EPS0 # F/m

def get_material_eps(material_name, w_cm_array):
    """Get 3x3 dielectric tensor for material at given wavenumbers"""
    if material_name == 'Air':
        eps = np.ones((len(w_cm_array), 3), dtype=complex)
    elif material_name == 'SiO2':
        # Simple SiO2 model - n=1.5 (can be improved with dispersion)
        n = 1.5
        eps = np.full((len(w_cm_array), 3), n**2, dtype=complex)
    else:
        # Try to get from builtin materials
        try:
            # Convert wavenumber to energy for builtin materials
            energy_eV = w_cm_array / 8065.54429  # cm⁻¹ to eV
            material = get_builtin_material(material_name, energy_eV=energy_eV[0])
            eps_tensor = material.get_dielectric_tensor(energy_eV[0])
            eps = np.full((len(w_cm_array), 3), eps_tensor, dtype=complex)
        except:
            raise ValueError(f"Unknown material: {material_name}")
    
    return eps.reshape(-1, 3, 1)

def build_layers(layer_info, w_cm_array):
    """Build layer structure for TMM calculation"""
    layers = []
    for material_name, thickness_nm in layer_info:
        eps = get_material_eps(material_name, w_cm_array)
        layers.append({
            'd': thickness_nm * 1e-9,  # Convert nm to m
            'eps': eps
        })
    return layers

def local_basis_rotation(kx, ky):
    """Calculate local coordinate system rotation matrix"""
    k_par = np.hypot(kx, ky)
    if k_par == 0:
        return np.eye(3)
    
    # s-polarization unit vector (perpendicular to k_par)
    s_hat = np.array([-ky, kx, 0]) / k_par
    z_hat = np.array([0, 0, 1])
    # p-polarization unit vector  
    p_hat = np.cross(z_hat, s_hat)
    
    return np.column_stack([s_hat, p_hat, z_hat])

def rotate_eps_full(eps_vec, kx, ky):
    """Rotate dielectric tensor to local coordinate system"""
    R = local_basis_rotation(kx, ky)
    eps_mat = np.diag(eps_vec)
    return R.T @ eps_mat @ R

def interface_jones(eps1, eps2, k_par, w):
    """Calculate Jones matrix for interface between two media"""
    # Calculate z-components of k-vectors for s and p polarizations
    kz1s = np.sqrt(eps1[1,1] * (w/c)**2 - k_par**2 + 0j)
    kz2s = np.sqrt(eps2[1,1] * (w/c)**2 - k_par**2 + 0j)
    kz1p = np.sqrt(eps1[0,0] * (w/c)**2 - (eps1[0,0]/eps1[2,2]) * k_par**2 + 0j)
    kz2p = np.sqrt(eps2[0,0] * (w/c)**2 - (eps2[0,0]/eps2[2,2]) * k_par**2 + 0j)
    
    # Calculate optical admittances
    Ys1 = kz1s / (mu0 * w)
    Yp1 = eps0 * eps1[2,2] * w / kz1p
    Ys2 = kz2s / (mu0 * w)  
    Yp2 = eps0 * eps2[2,2] * w / kz2p
    
    Y1 = np.diag([Ys1, Yp1])
    Y2 = np.diag([Ys2, Yp2])
    
    # Jones reflection matrix
    return np.linalg.solve(Y2 + Y1, Y2 - Y1)

def gamma_jones(rj, Gp, phi):
    """Propagate Jones matrix through layer with phase"""
    exp2 = np.exp(2j * phi)
    num = rj + Gp * exp2
    den = np.eye(2) + rj @ (Gp * exp2)
    return np.linalg.solve(den, num)

def calc_dispersion(layer_info, fixed_axis, fixed_val,
                   E_start, E_stop, q_start, q_stop,
                   a, b, delta_phi_deg):
    """
    Calculate 2D dispersion relations
    
    Parameters:
    - layer_info: List of (material_name, thickness_nm) tuples
    - fixed_axis: 'kx' or 'ky' - which axis to fix
    - fixed_val: Fixed value in μm⁻¹
    - E_start, E_stop: Energy range in eV
    - q_start, q_stop: q-vector range in μm⁻¹
    - a, b: Polarization amplitudes |Es|, |Ep|
    - delta_phi_deg: Phase difference in degrees
    
    Returns:
    - E_vals: Energy array (eV)
    - q_vals: q-vector array (μm⁻¹)
    - Rs, Rp: s and p reflectivities (2D arrays)
    - S0, S1, S2, S3: Stokes parameters (2D arrays)
    - phi, chi: Polarization angles (2D arrays)
    """
    # Resolution
    dE, dq = 0.01, 0.1
    E_vals = np.arange(E_start, E_stop + dE, dE)
    q_vals = np.arange(q_start, q_stop + dq, dq)
    nE, nQ = len(E_vals), len(q_vals)
    
    # Initialize result arrays
    Rs = np.zeros((nE, nQ))
    Rp = np.zeros((nE, nQ))
    S0 = np.zeros((nE, nQ))
    S1 = np.zeros((nE, nQ))
    S2 = np.zeros((nE, nQ))
    S3 = np.zeros((nE, nQ))
    phi = np.zeros((nE, nQ))
    chi = np.zeros((nE, nQ))
    
    print(f"Calculating dispersion: {nE} energies × {nQ} q-values = {nE*nQ} points")
    
    for i, E in enumerate(E_vals):
        if i % 50 == 0:
            print(f"Progress: {i}/{nE} energies ({100*i/nE:.1f}%)")
            
        # Convert energy to wavenumber and frequency
        w_spec = E * 8065.54429  # eV to cm⁻¹
        w = 2 * np.pi * c * w_spec * 100  # Angular frequency (rad/s)
        
        # Build layers for this energy
        layers = build_layers(layer_info, np.array([w_spec]))
        
        for j, q in enumerate(q_vals):
            # Set k-vectors based on fixed axis
            if fixed_axis == 'kx':
                kx = fixed_val * 1e6  # μm⁻¹ to m⁻¹
                ky = q * 1e6
            else:  # fixed_axis == 'ky'
                kx = q * 1e6
                ky = fixed_val * 1e6
                
            k_par = np.hypot(kx, ky)
            
            # Calculate Jones matrix for entire structure
            G = None
            for iface in reversed(range(len(layers) - 1)):
                # Get rotated dielectric tensors
                e1 = rotate_eps_full(layers[iface]['eps'][0, :, 0], kx, ky)
                e2 = rotate_eps_full(layers[iface + 1]['eps'][0, :, 0], kx, ky)
                
                # Interface Jones matrix
                rj = interface_jones(e1, e2, k_par, w)
                
                if iface == len(layers) - 2:  # Last interface
                    G = rj
                else:
                    # Propagation phase in layer
                    phi_prop = np.sqrt(e2[0,0] * (w/c)**2 - (e2[0,0]/e2[2,2]) * k_par**2 + 0j) * layers[iface + 1]['d']
                    G = gamma_jones(rj, G, phi_prop)
            
            # Input polarization in lab frame
            delta = np.deg2rad(delta_phi_deg)
            Ein_lab = np.array([-b * np.exp(1j * delta), a, 0])
            
            # Rotate to local frame
            R3 = local_basis_rotation(kx, ky)
            Ein_loc = R3.T @ Ein_lab
            
            # Apply Jones matrix
            if G is not None:
                Eout_loc = G @ Ein_loc[:2]
                
                # Rotate back to lab frame
                Eref_lab = R3 @ np.array([Eout_loc[0], Eout_loc[1], 0])
                Ex, Ey = Eref_lab[0], Eref_lab[1]
                
                # Calculate Stokes parameters
                S0[i,j] = abs(Ex)**2 + abs(Ey)**2
                if S0[i,j] > 1e-10:  # Avoid division by zero
                    S1[i,j] = (abs(Ex)**2 - abs(Ey)**2) / S0[i,j]
                    S2[i,j] = 2 * np.real(Ex * np.conj(Ey)) / S0[i,j]
                    S3[i,j] = 2 * np.imag(Ex * np.conj(Ey)) / S0[i,j]
                    
                    # Polarization angles
                    phi[i,j] = 0.5 * np.arctan2(S2[i,j], S1[i,j])
                    chi[i,j] = 0.5 * np.arcsin(np.clip(S3[i,j], -1, 1))
            
            # Reflectivities
            if G is not None:
                Rs[i,j] = abs(G[0,0])**2
                Rp[i,j] = abs(G[1,1])**2
    
    print("Dispersion calculation complete!")
    return E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi

def plot_dispersion_results(E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi, 
                          fixed_axis, fixed_val, save_figure=True):
    """Create 2D dispersion plots matching original Jones_disperion.py"""
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    
    # Plot parameters
    titles = ['Rs', 'Rp', 'S0', 'S1', 'S2', 'S3', 'phi', 'chi']
    data = [Rs, Rp, S0, S1, S2, S3, phi, chi]
    vmin_vals = [0.0, 0.0, -2.0, -1.0, -1.0, -1.0, -np.pi/2, -np.pi/4]
    vmax_vals = [1.0, 1.0,  2.0,  1.0,  1.0,  1.0,  np.pi/2,  np.pi/4]
    
    for idx, (ax, title, d) in enumerate(zip(axes.flatten(), titles, data)):
        im = ax.pcolormesh(q_vals, E_vals, d, shading='auto', cmap='viridis', 
                          vmin=vmin_vals[idx], vmax=vmax_vals[idx])
        ax.set_title(title)
        ax.set_xlabel('q (μm⁻¹)')
        ax.set_ylabel('E (eV)')
        fig.colorbar(im, ax=ax)
    
    # Add overall title
    fig.suptitle(f'2D Jones Dispersion Analysis (Fixed {fixed_axis} = {fixed_val} μm⁻¹)', 
                 fontsize=14, y=0.98)
    
    plt.tight_layout()
    
    if save_figure:
        filename = f'jones_dispersion_2d_{fixed_axis}_{fixed_val:.1f}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filename}")
    
    plt.show()

def main():
    """Main demonstration function"""
    print("=== 2D Jones Dispersion Analysis Demo ===")
    print("Reproducing original Jones_disperion.py functionality\n")
    
    # Layer structure (matching original default)
    layer_info = [
        ('Air', 1e6),      # Semi-infinite air (large thickness)
        ('SiO2', 500),     # 500 nm SiO2 layer
        ('Air', 1e6)       # Semi-infinite air substrate
    ]
    
    print("Layer structure:")
    for i, (material, thickness) in enumerate(layer_info):
        if thickness >= 1e6:
            print(f"  Layer {i+1}: {material} (semi-infinite)")
        else:
            print(f"  Layer {i+1}: {material} ({thickness} nm)")
    
    # Calculation parameters
    fixed_axis = 'kx'      # Fix kx, vary ky
    fixed_val = 0.0        # kx = 0 μm⁻¹ (normal incidence)
    E_start, E_stop = 1.0, 3.0    # Energy range (eV)
    q_start, q_stop = 0.0, 25.0   # q-vector range (μm⁻¹)
    
    # Polarization parameters
    a = 1.0          # |Es| amplitude
    b = 0.0          # |Ep| amplitude (s-polarized)
    delta_phi_deg = 0.0  # Phase difference
    
    print(f"\nCalculation parameters:")
    print(f"  Energy range: {E_start} - {E_stop} eV")
    print(f"  q-vector range: {q_start} - {q_stop} μm⁻¹")
    print(f"  Fixed axis: {fixed_axis} = {fixed_val} μm⁻¹")
    print(f"  Input polarization: |Es|={a}, |Ep|={b}, Δφ={delta_phi_deg}°")
    
    # Calculate dispersion
    print("\nStarting dispersion calculation...")
    results = calc_dispersion(layer_info, fixed_axis, fixed_val,
                             E_start, E_stop, q_start, q_stop,
                             a, b, delta_phi_deg)
    
    E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = results
    
    # Plot results
    print("\nCreating dispersion plots...")
    plot_dispersion_results(E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi,
                           fixed_axis, fixed_val, save_figure=True)
    
    print("\n=== Analysis Complete ===")
    print(f"Generated 2D dispersion plots with {len(E_vals)} × {len(q_vals)} data points")
    print("This reproduces the functionality of the original Jones_disperion.py file")

if __name__ == '__main__':
    main() 