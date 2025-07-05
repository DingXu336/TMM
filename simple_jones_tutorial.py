#!/usr/bin/env python3
"""
Simple Jones Calculator Tutorial

This tutorial shows you how to use the new unified Jones calculator
with step-by-step examples that actually work.

Author: Ding Xu (Physical Chemistry Researcher)
"""

import numpy as np
import matplotlib.pyplot as plt
from tmm.core.layer import Layer, create_substrate, create_superstrate, create_film
from tmm.materials.builtin_materials import get_air, get_sio2, get_au
from tmm.calculations.jones_3d_analyzer import Jones3DAnalyzer
from tmm.calculations.jones_unified_calculator import CalculationParams, GridResolution

def tutorial_step1_basic_structure():
    """Step 1: Create a basic layer structure"""
    print("="*60)
    print("STEP 1: Creating Basic Layer Structure")
    print("="*60)
    
    # Define energy range for materials (this is needed for material creation)
    energy_range = np.linspace(1.0, 3.0, 100)  # 1-3 eV
    
    # Create materials at the energy range
    air = get_air(energy_range)
    sio2 = get_sio2(energy_range, n=1.46)  # n=1.46 for SiO2
    
    print(f"✓ Created Air material: {air.name}")
    print(f"✓ Created SiO2 material: {sio2.name}")
    
    # Create layer structure
    layers = [
        create_superstrate(air, "Air (top)"),        # Semi-infinite air on top
        create_film(sio2, 500e-9, "SiO2 500nm"),    # 500 nm SiO2 layer  
        create_substrate(air, "Air (bottom)")        # Semi-infinite air substrate
    ]
    
    print(f"✓ Created layer structure with {len(layers)} layers:")
    for i, layer in enumerate(layers):
        thickness_str = "∞" if layer.is_semi_infinite else f"{layer.thickness_nm:.0f} nm"
        print(f"  {i+1}. {layer.name} ({thickness_str})")
    
    return layers

def tutorial_step2_quick_calculation():
    """Step 2: Run a quick coarse calculation"""
    print("\n" + "="*60)
    print("STEP 2: Quick Coarse Calculation (Fast Preview)")
    print("="*60)
    
    # Create layer structure
    layers = tutorial_step1_basic_structure()
    
    # Create the analyzer
    analyzer = Jones3DAnalyzer(layers)
    
    # Set up coarse parameters for quick calculation
    params = CalculationParams(
        energy_start_eV=1.5,
        energy_stop_eV=2.5,
        energy_points=11,        # Just 11 energy points
        kx_start_um=0.0,
        kx_stop_um=20.0,
        kx_points=21,           # 21 kx points
        ky_start_um=0.0,
        ky_stop_um=20.0,
        ky_points=21,           # 21 ky points
        es_amplitude=1.0,       # s-polarized light
        ep_amplitude=0.0,
        phase_diff_deg=0.0,
        resolution=GridResolution.COARSE
    )
    
    # Estimate memory usage
    memory_mb = analyzer.estimate_memory_usage(params)
    total_points = params.energy_points * params.kx_points * params.ky_points
    
    print(f"Grid size: {params.energy_points} × {params.kx_points} × {params.ky_points}")
    print(f"Total points: {total_points}")
    print(f"Estimated memory: {memory_mb:.1f} MB")
    print(f"Expected time: ~10-30 seconds")
    
    print("\nRunning calculation...")
    try:
        results = analyzer.calculate_3d_full(params)
        print(f"✓ Calculation completed in {results.computation_time:.1f} seconds!")
        return results, analyzer
    except Exception as e:
        print(f"✗ Calculation failed: {str(e)}")
        return None, None

def tutorial_step3_visualize_results(results, analyzer):
    """Step 3: Visualize the results"""
    print("\n" + "="*60)
    print("STEP 3: Visualizing Results")
    print("="*60)
    
    if results is None:
        print("No results to visualize")
        return
    
    # Plot results at middle energy
    mid_energy_idx = len(results.energy_eV) // 2
    energy_val = results.energy_eV[mid_energy_idx]
    
    print(f"Plotting results at energy: {energy_val:.2f} eV")
    
    # Create 2x2 plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot Rs (s-polarized reflectivity)
    im1 = axes[0, 0].pcolormesh(results.kx_um, results.ky_um, 
                               results.Rs[mid_energy_idx, :, :], 
                               shading='auto', cmap='viridis', vmin=0, vmax=1)
    axes[0, 0].set_title(f'Rs (s-reflectivity) at {energy_val:.2f} eV')
    axes[0, 0].set_xlabel('kx (μm⁻¹)')
    axes[0, 0].set_ylabel('ky (μm⁻¹)')
    plt.colorbar(im1, ax=axes[0, 0])
    
    # Plot Rp (p-polarized reflectivity)
    im2 = axes[0, 1].pcolormesh(results.kx_um, results.ky_um, 
                               results.Rp[mid_energy_idx, :, :], 
                               shading='auto', cmap='viridis', vmin=0, vmax=1)
    axes[0, 1].set_title(f'Rp (p-reflectivity) at {energy_val:.2f} eV')
    axes[0, 1].set_xlabel('kx (μm⁻¹)')
    axes[0, 1].set_ylabel('ky (μm⁻¹)')
    plt.colorbar(im2, ax=axes[0, 1])
    
    # Plot S0 (total intensity)
    im3 = axes[1, 0].pcolormesh(results.kx_um, results.ky_um, 
                               results.S0[mid_energy_idx, :, :], 
                               shading='auto', cmap='viridis')
    axes[1, 0].set_title(f'S0 (total intensity) at {energy_val:.2f} eV')
    axes[1, 0].set_xlabel('kx (μm⁻¹)')
    axes[1, 0].set_ylabel('ky (μm⁻¹)')
    plt.colorbar(im3, ax=axes[1, 0])
    
    # Plot S1 (linear polarization)
    im4 = axes[1, 1].pcolormesh(results.kx_um, results.ky_um, 
                               results.S1[mid_energy_idx, :, :], 
                               shading='auto', cmap='RdBu_r', vmin=-1, vmax=1)
    axes[1, 1].set_title(f'S1 (linear polarization) at {energy_val:.2f} eV')
    axes[1, 1].set_xlabel('kx (μm⁻¹)')
    axes[1, 1].set_ylabel('ky (μm⁻¹)')
    plt.colorbar(im4, ax=axes[1, 1])
    
    plt.suptitle('Basic SiO2 Layer Analysis (Tutorial Results)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('jones_tutorial_basic_results.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print("✓ Saved plot as: jones_tutorial_basic_results.png")

def tutorial_step4_extract_slices(results, analyzer):
    """Step 4: Extract 2D slices and 1D dispersion"""
    print("\n" + "="*60)
    print("STEP 4: Extracting 2D Slices and 1D Dispersion")
    print("="*60)
    
    if results is None:
        print("No results to extract from")
        return
    
    # Extract 2D slice at specific energy
    slice_energy = 2.0  # eV
    slice_2d = analyzer.get_2d_slice(results, slice_energy)
    
    print(f"✓ Extracted 2D k-space slice at {slice_2d['energy_eV']:.2f} eV")
    
    # Extract 1D dispersion along kx=0
    disp_1d = analyzer.get_1d_dispersion(results, 'kx', 0.0)
    
    print(f"✓ Extracted 1D dispersion along kx=0.0 μm⁻¹")
    
    # Plot the extracted data
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 2D slice
    im1 = ax1.pcolormesh(slice_2d['kx_um'], slice_2d['ky_um'], slice_2d['Rs'], 
                        shading='auto', cmap='viridis', vmin=0, vmax=1)
    ax1.set_title(f'2D k-space slice: Rs at {slice_2d["energy_eV"]:.2f} eV')
    ax1.set_xlabel('kx (μm⁻¹)')
    ax1.set_ylabel('ky (μm⁻¹)')
    plt.colorbar(im1, ax=ax1)
    
    # Plot 1D dispersion
    im2 = ax2.pcolormesh(disp_1d['k_um'], disp_1d['energy_eV'], disp_1d['Rs'], 
                        shading='auto', cmap='viridis', vmin=0, vmax=1)
    ax2.set_title('1D Dispersion: Rs vs Energy (kx=0)')
    ax2.set_xlabel('ky (μm⁻¹)')
    ax2.set_ylabel('Energy (eV)')
    plt.colorbar(im2, ax=ax2)
    
    plt.tight_layout()
    plt.savefig('jones_tutorial_extracted_data.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print("✓ Saved extraction plots as: jones_tutorial_extracted_data.png")

def tutorial_step5_plasmonic_example():
    """Step 5: More advanced example with plasmonic structure"""
    print("\n" + "="*60)
    print("STEP 5: Advanced Example - Plasmonic Structure")
    print("="*60)
    
    # Create plasmonic structure
    energy_range = np.linspace(0.5, 4.0, 100)  # Wider energy range for plasmonics
    
    air = get_air(energy_range)
    gold = get_au(energy_range)  # Gold with Drude model
    sio2 = get_sio2(energy_range, n=1.46)
    
    layers = [
        create_superstrate(air, "Air (top)"),
        create_film(gold, 50e-9, "Gold 50nm"),      # Thin gold layer
        create_substrate(sio2, "SiO2 (substrate)")
    ]
    
    print(f"✓ Created plasmonic structure:")
    for i, layer in enumerate(layers):
        thickness_str = "∞" if layer.is_semi_infinite else f"{layer.thickness_nm:.0f} nm"
        print(f"  {i+1}. {layer.name} ({thickness_str})")
    
    # Set up calculation for p-polarized light (better for plasmons)
    analyzer = Jones3DAnalyzer(layers)
    
    params = CalculationParams(
        energy_start_eV=1.0,
        energy_stop_eV=3.0,
        energy_points=16,        # Moderate resolution
        kx_start_um=0.0,
        kx_stop_um=25.0,
        kx_points=26,
        ky_start_um=0.0,
        ky_stop_um=25.0,
        ky_points=26,
        es_amplitude=0.0,       # No s-polarization
        ep_amplitude=1.0,       # p-polarized for plasmons
        phase_diff_deg=0.0,
        resolution=GridResolution.MEDIUM
    )
    
    memory_mb = analyzer.estimate_memory_usage(params)
    print(f"Memory estimate: {memory_mb:.1f} MB")
    
    print("\nRunning plasmonic calculation...")
    try:
        results = analyzer.calculate_3d_full(params)
        print(f"✓ Plasmonic calculation completed in {results.computation_time:.1f} seconds!")
        
        # Plot Rp results at middle energy
        mid_energy_idx = len(results.energy_eV) // 2
        energy_val = results.energy_eV[mid_energy_idx]
        
        plt.figure(figsize=(8, 6))
        plt.pcolormesh(results.kx_um, results.ky_um, 
                      results.Rp[mid_energy_idx, :, :], 
                      shading='auto', cmap='viridis', vmin=0, vmax=1)
        plt.colorbar(label='Rp (p-reflectivity)')
        plt.title(f'Plasmonic Structure: Rp at {energy_val:.2f} eV')
        plt.xlabel('kx (μm⁻¹)')
        plt.ylabel('ky (μm⁻¹)')
        plt.savefig('jones_tutorial_plasmonic.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        print("✓ Saved plasmonic plot as: jones_tutorial_plasmonic.png")
        
    except Exception as e:
        print(f"✗ Plasmonic calculation failed: {str(e)}")

def main():
    """Main tutorial function"""
    print("="*60)
    print("JONES CALCULATOR TUTORIAL")
    print("="*60)
    print("This tutorial will teach you how to use the unified Jones calculator.")
    print("We'll go through 5 steps from basic to advanced usage.")
    print("="*60)
    
    try:
        # Step 1: Create basic structure (this just shows how, doesn't run calculation)
        tutorial_step1_basic_structure()
        
        # Step 2: Run quick calculation
        results, analyzer = tutorial_step2_quick_calculation()
        
        # Step 3: Visualize results
        tutorial_step3_visualize_results(results, analyzer)
        
        # Step 4: Extract data
        tutorial_step4_extract_slices(results, analyzer)
        
        # Step 5: Advanced example
        tutorial_step5_plasmonic_example()
        
        print("\n" + "="*60)
        print("✓ TUTORIAL COMPLETED SUCCESSFULLY!")
        print("="*60)
        print("Generated files:")
        print("• jones_tutorial_basic_results.png")
        print("• jones_tutorial_extracted_data.png") 
        print("• jones_tutorial_plasmonic.png")
        print("\nNext steps:")
        print("• Try the GUI: python launch_jones_gui.py")
        print("• Modify parameters to explore different structures")
        print("• Use finer grids for higher quality results")
        print("="*60)
        
    except Exception as e:
        print(f"\n✗ Tutorial failed: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 