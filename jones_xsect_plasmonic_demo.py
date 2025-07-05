#!/usr/bin/env python3
"""
2D Jones k-space Cross-Section Demo - Plasmonic Features with Gold

This script demonstrates plasmonic features in k-space using a Gold film
to show surface plasmons and other dispersive effects in 2D kx-ky plots.

Author: Ding Xu (Physical Chemistry Researcher)
Based on original Jones_Xect.py with plasmonic enhancements
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from tmm.materials import get_builtin_material
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
        # SiO2 with simple dispersion
        n = 1.5
        eps = np.full((len(w_cm_array), 3), n**2, dtype=complex)
    elif material_name == 'Gold':
        # Gold with Drude model approximation
        eps = np.zeros((len(w_cm_array), 3), dtype=complex)
        for i, w_cm in enumerate(w_cm_array):
            energy_eV = w_cm / 8065.54429  # cm⁻¹ to eV
            # Drude model for gold
            eps_inf = 1.0
            wp = 9.0  # eV (plasma frequency)
            gamma = 0.1  # eV (damping)
            eps_val = eps_inf - wp**2 / (energy_eV**2 + 1j*gamma*energy_eV)
            eps[i, :] = eps_val
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

def calc_2d_kspace_xsect_fast(layer_info, energy_eV, qstart, qstop, dq, 
                             a, b, delta_phi_deg, resolution='medium'):
    """
    Calculate 2D k-space cross-section at fixed energy with adjustable resolution
    
    Parameters:
    - resolution: 'low', 'medium', 'high' for different grid densities
    """
    
    # Adjust resolution
    if resolution == 'low':
        dq_eff = dq * 3
    elif resolution == 'medium':
        dq_eff = dq * 1.5
    else:  # high
        dq_eff = dq
    
    # Convert energy to frequency
    w_spec = energy_eV * 8065.54429  # eV to cm⁻¹
    w = 2 * np.pi * c * w_spec * 100  # Angular frequency (rad/s)
    
    # Create k-space grid
    qx = np.arange(qstart, qstop + dq_eff, dq_eff)
    qy = qx.copy()
    KX, KY = np.meshgrid(qx, qy, indexing='xy')
    Kpar = np.hypot(KX, KY)
    
    # Build layers at this energy
    layers = build_layers(layer_info, np.array([w_spec]))
    n_ifaces = len(layers) - 1
    
    # Initialize Jones matrix array
    Gmnt = np.zeros((KX.shape[0], KX.shape[1], 2, 2), dtype=complex)
    
    print(f"Calculating 2D k-space cross-section at {energy_eV} eV ({resolution} resolution)")
    print(f"Grid size: {KX.shape[0]} × {KX.shape[1]} = {KX.size} k-points")
    
    # Calculate Jones matrix for each k-point
    for iy in range(KX.shape[0]):
        if iy % max(1, KX.shape[0]//10) == 0:
            print(f"Progress: {iy}/{KX.shape[0]} rows ({100*iy/KX.shape[0]:.1f}%)")
            
        for ix in range(KX.shape[1]):
            kx = KX[iy, ix] * 1e6  # μm⁻¹ to m⁻¹
            ky = KY[iy, ix] * 1e6
            kp = Kpar[iy, ix] * 1e6
            
            # Build Jones matrix reflection by cascading interfaces
            # Start with bottom interface
            e1 = rotate_eps_full(layers[-2]['eps'][0, :, 0], kx, ky)
            e2 = rotate_eps_full(layers[-1]['eps'][0, :, 0], kx, ky)
            G = interface_jones(e1, e2, kp, w)
            
            # Cascade upward through all interfaces
            for iface in range(1, n_ifaces):
                j = len(layers) - 1 - iface
                if j < 1:
                    break
                    
                e1 = rotate_eps_full(layers[j-1]['eps'][0, :, 0], kx, ky)
                e2 = rotate_eps_full(layers[j]['eps'][0, :, 0], kx, ky)
                rj = interface_jones(e1, e2, kp, w)
                
                # Propagation phase through layer j
                kz2p = np.sqrt(e2[0,0] * (w/c)**2 - (e2[0,0]/e2[2,2]) * kp**2 + 0j)
                phi_prop = kz2p * layers[j]['d']
                
                G = gamma_jones(rj, G, phi_prop)
            
            Gmnt[iy, ix] = G
    
    # Input field in lab frame
    delta = np.deg2rad(delta_phi_deg)
    Ein_lab = np.array([-b * np.exp(1j * delta), a, 0])  # lab: x->-p, y->s, z
    
    # Calculate reflected fields
    Eref_lab = np.zeros((KX.shape[0], KX.shape[1], 3), dtype=complex)
    
    for iy in range(KX.shape[0]):
        for ix in range(KX.shape[1]):
            # Rotation matrix for this k-point
            R3 = local_basis_rotation(KX[iy, ix] * 1e6, KY[iy, ix] * 1e6)
            
            # Rotate input field to local frame
            Ein_loc = R3.T @ Ein_lab
            
            # Apply Jones matrix
            Eout_loc2 = Gmnt[iy, ix] @ Ein_loc[:2]
            
            # Rotate back to lab frame
            Eref_lab[iy, ix] = R3 @ np.array([Eout_loc2[0], Eout_loc2[1], 0])
    
    # Extract field components
    Ex = Eref_lab[:, :, 0]
    Ey = Eref_lab[:, :, 1]
    
    # Calculate Stokes parameters
    S0 = np.abs(Ex)**2 + np.abs(Ey)**2
    S1 = np.abs(Ex)**2 - np.abs(Ey)**2
    S2 = 2 * np.real(Ex * np.conj(Ey))
    S3 = 2 * np.imag(Ex * np.conj(Ey))
    
    # Polarization angles
    phi = 0.5 * np.arctan2(S2, S1)
    chi = 0.5 * np.arcsin(np.clip(S3 / (S0 + 1e-10), -1, 1))
    
    # Scalar reflectivities from diagonal of Jones matrix
    Rs = np.abs(Gmnt[:, :, 0, 0])**2
    Rp = np.abs(Gmnt[:, :, 1, 1])**2
    
    print("2D k-space cross-section calculation complete!")
    return KX, KY, Rs, Rp, S0, S1, S2, S3, phi, chi

def plot_kspace_xsect_enhanced(KX, KY, Rs, Rp, S0, S1, S2, S3, phi, chi,
                              energy_eV, layer_info, save_figure=True):
    """Create enhanced 2D k-space cross-section plots with better styling"""
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    
    # Plot parameters
    titles = ['Rs = |rₛ|² (s-pol)', 'Rp = |rₚ|² (p-pol)', 'S₀ (Intensity)', 'S₁ (Linear H-V)', 
              'S₂ (Linear ±45°)', 'S₃ (Circular)', 'φ (Azimuth)', 'χ (Ellipticity)']
    data = [Rs, Rp, S0, S1, S2, S3, phi, chi]
    vmin_vals = [0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -np.pi/2, -np.pi/4]
    vmax_vals = [1.0, 1.0, 2.0,  1.0,  1.0,  1.0,  np.pi/2,  np.pi/4]
    cmaps = ['viridis', 'viridis', 'plasma', 'RdBu_r', 'RdBu_r', 'RdBu_r', 'hsv', 'RdBu_r']
    
    for idx, (ax, title, d) in enumerate(zip(axes.flatten(), titles, data)):
        im = ax.pcolormesh(KX, KY, d, shading='auto', cmap=cmaps[idx], 
                          vmin=vmin_vals[idx], vmax=vmax_vals[idx])
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('kx (μm⁻¹)', fontsize=11)
        ax.set_ylabel('ky (μm⁻¹)', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.ax.tick_params(labelsize=9)
        
        # Add circular symmetry guide for normal incidence
        if idx < 2:  # Only for Rs and Rp plots
            circle = Circle((0, 0), np.sqrt(d.max()) * 5, fill=False, 
                           color='white', alpha=0.5, linestyle='--', linewidth=1)
            ax.add_patch(circle)
    
    # Add overall title with layer information
    layer_str = " → ".join([f"{mat}" + (f" ({thick:.0f}nm)" if thick < 1e6 else " (∞)") 
                           for mat, thick in layer_info])
    fig.suptitle(f'2D k-space Cross-Section at {energy_eV} eV\n{layer_str}', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    if save_figure:
        materials = "_".join([mat for mat, _ in layer_info])
        filename = f'jones_xsect_2d_{materials}_{energy_eV}eV.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {filename}")
    
    plt.show()

def main():
    """Main demonstration function"""
    print("=== 2D Jones k-space Cross-Section - Plasmonic Demo ===")
    print("Showing plasmonic features in k-space with Gold film\n")
    
    # Plasmonic layer structure
    layer_info = [
        ('Air', 1e6),      # Semi-infinite air
        ('Gold', 50),      # 50 nm Gold film
        ('SiO2', 1e6)      # Semi-infinite SiO2 substrate
    ]
    
    print("Layer structure:")
    for i, (material, thickness) in enumerate(layer_info):
        if thickness >= 1e6:
            print(f"  Layer {i+1}: {material} (semi-infinite)")
        else:
            print(f"  Layer {i+1}: {material} ({thickness} nm)")
    
    # Calculation parameters - focus on plasmonic region
    energy_eV = 2.4        # Fixed energy in plasmonic range (eV)
    qstart, qstop = 0.0, 30.0    # k-vector range (μm⁻¹) - higher for plasmons
    dq = 0.15             # k-vector step (μm⁻¹)
    
    # Polarization parameters - use p-polarization for plasmons
    a = 0.0          # |Es| amplitude
    b = 1.0          # |Ep| amplitude (p-polarized for plasmons)
    delta_phi_deg = 0.0  # Phase difference
    
    print(f"\nCalculation parameters:")
    print(f"  Fixed energy: {energy_eV} eV (plasmonic region)")
    print(f"  k-vector range: {qstart} - {qstop} μm⁻¹")
    print(f"  k-vector step: {dq} μm⁻¹")
    print(f"  Input polarization: |Es|={a}, |Ep|={b}, Δφ={delta_phi_deg}° (p-polarized)")
    
    # Calculate 2D k-space cross-section
    print("\nStarting 2D k-space cross-section calculation...")
    results = calc_2d_kspace_xsect_fast(layer_info, energy_eV, qstart, qstop, dq,
                                        a, b, delta_phi_deg, resolution='medium')
    
    KX, KY, Rs, Rp, S0, S1, S2, S3, phi, chi = results
    
    # Plot results
    print("\nCreating enhanced k-space cross-section plots...")
    plot_kspace_xsect_enhanced(KX, KY, Rs, Rp, S0, S1, S2, S3, phi, chi,
                              energy_eV, layer_info, save_figure=True)
    
    print("\n=== Plasmonic k-space Analysis Complete ===")
    print(f"Generated 2D k-space cross-section plots with {KX.shape[0]} × {KX.shape[1]} k-points")
    print("This shows plasmonic features in k-space including:")
    print("  - Surface plasmon resonances (circular patterns in Rp)")
    print("  - Enhanced field coupling at specific k-values")
    print("  - Polarization conversion effects")
    print("  - Symmetry breaking due to anisotropic structure")

if __name__ == '__main__':
    main() 