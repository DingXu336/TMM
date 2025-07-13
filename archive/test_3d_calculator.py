#!/usr/bin/env python3
"""
Test script for the 3D TMM Calculator
Demonstrates basic usage with coarse grid for testing
"""

import numpy as np
import matplotlib.pyplot as plt
from Jones_3D_Calculator import calc_3d_reflectance, TMM3DGui
import time

def test_3d_calculation():
    """Test the 3D calculation function with a simple structure"""
    print("Testing 3D TMM calculation...")
    
    # Define a simple three-layer structure
    layer_info = [
        ('Air', 1e6),    # Semi-infinite substrate
        ('SiO2', 500),   # 500 nm SiO2 layer
        ('Air', 1e6)     # Semi-infinite superstrate
    ]
    
    # Coarse grid parameters for testing
    E_start, E_stop, dE = 1.0, 2.0, 0.5  # 3 energy points
    kx_start, kx_stop, dkx = -5, 5, 5     # 3 kx points
    ky_start, ky_stop, dky = -5, 5, 5     # 3 ky points
    
    # Input polarization (s-polarized)
    a, b = 1.0, 0.0  # |Es|=1, |Ep|=0
    delta_phi_deg = 0.0
    
    print(f"Grid size: {int((E_stop-E_start)/dE)+1} x {int((kx_stop-kx_start)/dkx)+1} x {int((ky_stop-ky_start)/dky)+1}")
    print("Calculating...")
    
    start_time = time.time()
    
    # Progress callback for testing
    def progress_callback(progress):
        print(f"Progress: {progress*100:.1f}%")
    
    # Perform calculation
    results = calc_3d_reflectance(
        layer_info, E_start, E_stop, dE,
        kx_start, kx_stop, dkx, ky_start, ky_stop, dky,
        a, b, delta_phi_deg, progress_callback
    )
    
    end_time = time.time()
    
    E_vals, kx_vals, ky_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = results
    
    print(f"Calculation completed in {end_time - start_time:.2f} seconds")
    print(f"Results shape: {Rs.shape}")
    print(f"Energy range: {E_vals[0]:.2f} - {E_vals[-1]:.2f} eV")
    print(f"kx range: {kx_vals[0]:.2f} - {kx_vals[-1]:.2f} μm⁻¹")
    print(f"ky range: {ky_vals[0]:.2f} - {ky_vals[-1]:.2f} μm⁻¹")
    print(f"Rs range: {Rs.min():.3f} - {Rs.max():.3f}")
    print(f"Rp range: {Rp.min():.3f} - {Rp.max():.3f}")
    
    # Simple visualization
    plt.figure(figsize=(12, 8))
    
    # Plot Rs for middle energy
    E_mid = len(E_vals) // 2
    plt.subplot(2, 2, 1)
    KX, KY = np.meshgrid(kx_vals, ky_vals, indexing='ij')
    plt.pcolormesh(KX, KY, Rs[E_mid, :, :], shading='auto', cmap='viridis')
    plt.colorbar(label='Rs')
    plt.xlabel('kx (μm⁻¹)')
    plt.ylabel('ky (μm⁻¹)')
    plt.title(f'Rs at E = {E_vals[E_mid]:.2f} eV')
    
    # Plot Rp for middle energy
    plt.subplot(2, 2, 2)
    plt.pcolormesh(KX, KY, Rp[E_mid, :, :], shading='auto', cmap='viridis')
    plt.colorbar(label='Rp')
    plt.xlabel('kx (μm⁻¹)')
    plt.ylabel('ky (μm⁻¹)')
    plt.title(f'Rp at E = {E_vals[E_mid]:.2f} eV')
    
    # Plot S1 for middle energy
    plt.subplot(2, 2, 3)
    plt.pcolormesh(KX, KY, S1[E_mid, :, :], shading='auto', cmap='RdBu', vmin=-1, vmax=1)
    plt.colorbar(label='S1')
    plt.xlabel('kx (μm⁻¹)')
    plt.ylabel('ky (μm⁻¹)')
    plt.title(f'S1 at E = {E_vals[E_mid]:.2f} eV')
    
    # Plot energy dispersion at center of k-space
    plt.subplot(2, 2, 4)
    kx_mid, ky_mid = len(kx_vals) // 2, len(ky_vals) // 2
    plt.plot(E_vals, Rs[:, kx_mid, ky_mid], 'b-', label='Rs')
    plt.plot(E_vals, Rp[:, kx_mid, ky_mid], 'r-', label='Rp')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Reflectance')
    plt.title(f'Energy dispersion at kx={kx_vals[kx_mid]:.1f}, ky={ky_vals[ky_mid]:.1f}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return results

def test_gui():
    """Test the GUI application"""
    print("Starting GUI application...")
    print("Use the 'Coarse (Test)' preset for quick testing")
    print("Use 'Medium' or 'Fine (Final)' for production calculations")
    
    app = TMM3DGui()
    app.run()

if __name__ == '__main__':
    import sys
    
    print("3D TMM Calculator Test Script")
    print("=" * 40)
    print("1. Command-line calculation test")
    print("2. GUI application test")
    print("=" * 40)
    
    choice = input("Enter choice (1 or 2): ").strip()
    
    if choice == '1':
        test_3d_calculation()
    elif choice == '2':
        test_gui()
    else:
        print("Invalid choice. Running calculation test by default...")
        test_3d_calculation() 