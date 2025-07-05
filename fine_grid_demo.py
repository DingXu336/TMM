#!/usr/bin/env python3
"""
Fine Grid Demo: Air(100μm)/SiO2(500nm)/Air(100μm)

This demo shows high-resolution 3D E-kx-ky analysis of a simple
dielectric layer structure using fine grid settings.

Author: Ding Xu (Physical Chemistry Researcher)
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import traceback

from tmm.core.layer import Layer, create_superstrate, create_film, create_substrate
from tmm.materials.builtin_materials import get_builtin_material
from tmm.calculations.jones_3d_analyzer import Jones3DAnalyzer
from tmm.calculations.jones_unified_calculator import CalculationParams, GridResolution

def create_air_sio2_air_structure():
    """Create Air(100μm)/SiO2(500nm)/Air(100μm) structure"""
    
    # Create materials
    air = get_builtin_material("Air")
    sio2 = get_builtin_material("SiO2")
    
    # Create layer structure
    layers = [
        create_superstrate(air),                      # Semi-infinite air superstrate
        create_film(sio2, thickness_nm=500),          # 500 nm SiO2 film
        create_substrate(air)                         # Semi-infinite air substrate
    ]
    
    return layers

def main():
    """Run fine grid demo"""
    print("="*80)
    print("FINE GRID DEMO: Air(100μm)/SiO2(500nm)/Air(100μm)")
    print("="*80)
    print()
    
    try:
        # Create layer structure
        print("Creating layer structure...")
        layers = create_air_sio2_air_structure()
        
        # Display structure
        print("Layer structure:")
        for i, layer in enumerate(layers):
            if layer.is_semi_infinite:
                print(f"  Layer {i+1}: {layer.material.name} - semi-infinite")
            else:
                print(f"  Layer {i+1}: {layer.material.name} - {layer.thickness_nm:.0f} nm")
        print()
        
        # Create fine grid parameters
        print("Setting up fine grid parameters...")
        params = CalculationParams(
            # Energy range
            energy_start_eV=1.0,
            energy_stop_eV=3.0,
            energy_points=61,  # Fine energy resolution
            
            # k-space range
            kx_start_um=0.0,
            kx_stop_um=25.0,
            kx_points=51,      # Fine kx resolution
            
            ky_start_um=0.0,
            ky_stop_um=25.0,
            ky_points=51,      # Fine ky resolution
            
            # Polarization (s-polarized)
            es_amplitude=1.0,
            ep_amplitude=0.0,
            phase_diff_deg=0.0,
            
            # Fine grid resolution
            resolution=GridResolution.FINE,
            max_memory_gb=4.0
        )
        
        # Display parameters
        print("Calculation parameters:")
        print(f"  Energy: {params.energy_start_eV:.1f} - {params.energy_stop_eV:.1f} eV ({params.energy_points} points)")
        print(f"  kx: {params.kx_start_um:.1f} - {params.kx_stop_um:.1f} μm⁻¹ ({params.kx_points} points)")
        print(f"  ky: {params.ky_start_um:.1f} - {params.ky_stop_um:.1f} μm⁻¹ ({params.ky_points} points)")
        print(f"  Total grid points: {params.energy_points * params.kx_points * params.ky_points:,}")
        print(f"  Grid resolution: {params.resolution.value}")
        print(f"  Polarization: |Es|={params.es_amplitude:.1f}, |Ep|={params.ep_amplitude:.1f}")
        print()
        
        # Estimate memory usage
        total_points = params.energy_points * params.kx_points * params.ky_points
        memory_gb = total_points * 8 * 8 / (1024**3)  # 8 quantities × 8 bytes per complex number
        print(f"Estimated memory usage: {memory_gb:.2f} GB")
        print()
        
        # Create analyzer
        print("Creating 3D analyzer...")
        analyzer = Jones3DAnalyzer(layers)
        
        # Run calculation with timing
        print("Starting fine grid calculation...")
        print("This may take several minutes due to fine resolution...")
        start_time = time.time()
        
        # Progress callback
        def progress_callback(completed, total, eta_seconds=None):
            percent = (completed / total) * 100
            if eta_seconds is not None:
                eta_min = eta_seconds / 60
                print(f"  Progress: {completed}/{total} ({percent:.1f}%) - ETA: {eta_min:.1f} minutes")
            else:
                print(f"  Progress: {completed}/{total} ({percent:.1f}%)")
        
        # Calculate
        results = analyzer.calculate_3d_full(params)
        
        calc_time = time.time() - start_time
        print(f"Calculation completed in {calc_time:.1f} seconds ({calc_time/60:.1f} minutes)")
        print()
        
        # Display results summary
        print("Results summary:")
        print(f"  Data shape: {results.Rs.shape}")
        print(f"  Available quantities: ['Rs', 'Rp', 'S0', 'S1', 'S2', 'S3', 'phi', 'chi']")
        print(f"  Energy range: {results.energy_eV[0]:.2f} - {results.energy_eV[-1]:.2f} eV")
        print(f"  kx range: {results.kx_um[0]:.2f} - {results.kx_um[-1]:.2f} μm⁻¹")
        print(f"  ky range: {results.ky_um[0]:.2f} - {results.ky_um[-1]:.2f} μm⁻¹")
        print()
        
        # Create visualization
        print("Creating visualization...")
        
        # Create figure with 2x4 subplots
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        fig.suptitle('Fine Grid Demo: Air(100μm)/SiO2(500nm)/Air(100μm)\n'
                    f'E-kx slice at ky=0 - {params.energy_points}×{params.kx_points} points', 
                    fontsize=16, fontweight='bold')
        
        # Extract 2D slice at ky=0
        ky_center_idx = params.ky_points // 2
        
        # Plot quantities
        quantities = ['Rs', 'Rp', 'S0', 'S1', 'S2', 'S3', 'phi', 'chi']
        titles = ['Rs (s-reflectivity)', 'Rp (p-reflectivity)', 'S0 (total intensity)', 'S1 (linear H-V)',
                 'S2 (linear ±45°)', 'S3 (circular)', 'φ (azimuth)', 'χ (ellipticity)']
        
        for i, (qty, title) in enumerate(zip(quantities, titles)):
            ax = axes[i//4, i%4]
            
            # Get 2D slice
            if hasattr(results, qty):
                data_slice = getattr(results, qty)[:, :, ky_center_idx]
                
                # Plot with appropriate colormap
                if qty in ['Rs', 'Rp', 'S0']:
                    im = ax.imshow(data_slice, aspect='auto', origin='lower', 
                                 extent=[results.kx_um[0], results.kx_um[-1], 
                                        results.energy_eV[0], results.energy_eV[-1]],
                                 cmap='hot', vmin=0, vmax=1)
                else:
                    im = ax.imshow(data_slice, aspect='auto', origin='lower',
                                 extent=[results.kx_um[0], results.kx_um[-1], 
                                        results.energy_eV[0], results.energy_eV[-1]],
                                 cmap='RdBu_r', vmin=-1, vmax=1)
                
                plt.colorbar(im, ax=ax, shrink=0.8)
                ax.set_title(title, fontsize=12)
                ax.set_xlabel('kx (μm⁻¹)')
                ax.set_ylabel('Energy (eV)')
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        filename = f'fine_grid_air_sio2_air_E{params.energy_points}_kx{params.kx_points}_ky{params.ky_points}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as: {filename}")
        
        # Show plot
        plt.show()
        
        # Extract some analysis
        print("\nAnalysis:")
        print(f"  Max Rs: {np.max(results.Rs):.4f}")
        print(f"  Min Rs: {np.min(results.Rs):.4f}")
        print(f"  Max Rp: {np.max(results.Rp):.4f}")
        print(f"  Min Rp: {np.min(results.Rp):.4f}")
        
        # Save data
        save_choice = input("\nSave calculation data? (y/n): ").lower().strip()
        if save_choice == 'y':
            import pickle
            data_filename = f'fine_grid_air_sio2_air_data.pkl'
            with open(data_filename, 'wb') as f:
                pickle.dump({
                    'results': results,
                    'params': params,
                    'layers': layers,
                    'calculation_time': calc_time
                }, f)
            print(f"Data saved as: {data_filename}")
        
        print("\nFine grid demo completed successfully!")
        
    except Exception as e:
        print(f"Error during calculation: {str(e)}")
        print("Traceback:")
        traceback.print_exc()
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    if success:
        print("\n" + "="*80)
        print("FINE GRID DEMO COMPLETED")
        print("="*80)
    else:
        print("\n" + "="*80)
        print("FINE GRID DEMO FAILED")
        print("="*80) 