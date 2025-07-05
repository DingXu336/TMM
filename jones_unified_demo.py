#!/usr/bin/env python3
"""
Unified Jones Calculator Demo

This script demonstrates the unified Jones calculator capabilities:
- 3D E-kx-ky analysis with adaptive grid resolution
- Memory-efficient calculations with coarse/fine grid options
- Enhanced visualization and data export

Features demonstrated:
1. Coarse grid for fast preview
2. Fine grid for detailed analysis
3. Memory estimation and management
4. Flexible output options (3D, 2D slices, 1D dispersion)

Author: Ding Xu (Physical Chemistry Researcher)
"""

import numpy as np
import matplotlib.pyplot as plt
from tmm.core.layer import create_substrate, create_superstrate, create_film
from tmm.materials.builtin_materials import get_builtin_material
from tmm.calculations.jones_unified_calculator import (
    JonesUnifiedCalculator, CalculationParams, GridResolution
)
from tmm.calculations.jones_3d_analyzer import Jones3DAnalyzer
import time

def create_basic_structure():
    """Create a basic SiO2 layer structure"""
    print("Creating basic SiO2 structure...")
    
    # Get materials
    air = get_builtin_material("Air")
    sio2 = get_builtin_material("SiO2")
    
    # Create layer structure
    layers = [
        create_superstrate(air, "Air (superstrate)"),
        create_film(sio2, 500e-9, "SiO2 500nm"),  # 500 nm in meters
        create_substrate(air, "Air (substrate)")
    ]
    
    return layers

def create_plasmonic_structure():
    """Create a plasmonic structure with gold layer"""
    print("Creating plasmonic structure...")
    
    # Get materials
    air = get_builtin_material("Air")
    gold = get_builtin_material("Gold")
    sio2 = get_builtin_material("SiO2")
    
    # Create layer structure
    layers = [
        create_superstrate(air, "Air (superstrate)"),
        create_film(gold, 50e-9, "Gold 50nm"),  # 50 nm in meters
        create_substrate(sio2, "SiO2 (substrate)")
    ]
    
    return layers

def demo_coarse_grid_preview():
    """Demo 1: Coarse grid for fast preview"""
    print("\n" + "="*60)
    print("DEMO 1: Coarse Grid Preview (Fast)")
    print("="*60)
    
    # Create structure
    layers = create_basic_structure()
    
    # Create analyzer
    analyzer = Jones3DAnalyzer(layers)
    
    # Set up coarse grid parameters
    params = CalculationParams(
        energy_start_eV=1.0,
        energy_stop_eV=3.0,
        energy_points=21,  # Coarse energy grid
        kx_start_um=0.0,
        kx_stop_um=25.0,
        kx_points=26,     # Coarse kx grid
        ky_start_um=0.0,
        ky_stop_um=25.0,
        ky_points=26,     # Coarse ky grid
        es_amplitude=1.0,
        ep_amplitude=0.0,
        phase_diff_deg=0.0,
        resolution=GridResolution.COARSE
    )
    
    # Estimate memory usage
    memory_mb = analyzer.estimate_memory_usage(params)
    print(f"Estimated memory usage: {memory_mb:.1f} MB")
    print(f"Grid size: {params.energy_points} × {params.kx_points} × {params.ky_points}")
    print(f"Total points: {params.energy_points * params.kx_points * params.ky_points}")
    
    # Run calculation
    print("\nRunning coarse grid calculation...")
    start_time = time.time()
    
    try:
        results = analyzer.calculate_3d_full(params)
        
        print(f"✓ Coarse calculation completed in {results.computation_time:.2f} seconds")
        print(f"✓ Memory usage: {results.memory_usage_mb:.1f} MB")
        
        # Plot preview
        plot_preview_results(results, "Coarse Grid Preview")
        
        return results
    
    except Exception as e:
        print(f"✗ Coarse calculation failed: {str(e)}")
        return None

def demo_fine_grid_analysis():
    """Demo 2: Fine grid for detailed analysis"""
    print("\n" + "="*60)
    print("DEMO 2: Fine Grid Analysis (Detailed)")
    print("="*60)
    
    # Create plasmonic structure for more interesting results
    layers = create_plasmonic_structure()
    
    # Create analyzer
    analyzer = Jones3DAnalyzer(layers)
    
    # Set up fine grid parameters
    params = CalculationParams(
        energy_start_eV=0.5,
        energy_stop_eV=4.0,
        energy_points=36,  # Fine energy grid
        kx_start_um=0.0,
        kx_stop_um=30.0,
        kx_points=31,     # Fine kx grid
        ky_start_um=0.0,
        ky_stop_um=30.0,
        ky_points=31,     # Fine ky grid
        es_amplitude=0.0,
        ep_amplitude=1.0,  # p-polarized for plasmons
        phase_diff_deg=0.0,
        resolution=GridResolution.FINE
    )
    
    # Estimate memory usage
    memory_mb = analyzer.estimate_memory_usage(params)
    print(f"Estimated memory usage: {memory_mb:.1f} MB")
    print(f"Grid size: {params.energy_points} × {params.kx_points} × {params.ky_points}")
    print(f"Total points: {params.energy_points * params.kx_points * params.ky_points}")
    
    # Check if calculation is too large
    if memory_mb > 500:  # 500 MB limit for demo
        print(f"⚠ Calculation too large for demo (>{500} MB). Using medium grid instead.")
        params.resolution = GridResolution.MEDIUM
        params.energy_points = 18
        params.kx_points = 21
        params.ky_points = 21
        memory_mb = analyzer.estimate_memory_usage(params)
        print(f"Adjusted memory usage: {memory_mb:.1f} MB")
    
    # Run calculation
    print("\nRunning fine grid calculation...")
    
    try:
        results = analyzer.calculate_3d_full(params)
        
        print(f"✓ Fine calculation completed in {results.computation_time:.2f} seconds")
        print(f"✓ Memory usage: {results.memory_usage_mb:.1f} MB")
        
        # Plot detailed results
        plot_detailed_results(results, "Fine Grid Analysis - Plasmonic Structure")
        
        return results
    
    except Exception as e:
        print(f"✗ Fine calculation failed: {str(e)}")
        return None

def demo_2d_slice_extraction():
    """Demo 3: 2D slice extraction from 3D results"""
    print("\n" + "="*60)
    print("DEMO 3: 2D Slice Extraction")
    print("="*60)
    
    # Create structure
    layers = create_basic_structure()
    analyzer = Jones3DAnalyzer(layers)
    
    # Quick 3D calculation for slicing
    params = CalculationParams(
        energy_start_eV=1.0,
        energy_stop_eV=3.0,
        energy_points=21,
        kx_start_um=0.0,
        kx_stop_um=25.0,
        kx_points=26,
        ky_start_um=0.0,
        ky_stop_um=25.0,
        ky_points=26,
        es_amplitude=1.0,
        ep_amplitude=0.0,
        phase_diff_deg=0.0,
        resolution=GridResolution.MEDIUM
    )
    
    print("Running 3D calculation for slice extraction...")
    try:
        results = analyzer.calculate_3d_full(params)
        
        # Extract 2D slices at different energies
        energies_to_slice = [1.5, 2.0, 2.5]
        
        plt.figure(figsize=(15, 5))
        
        for i, energy in enumerate(energies_to_slice):
            slice_data = analyzer.get_2d_slice(results, energy)
            
            plt.subplot(1, 3, i+1)
            plt.pcolormesh(slice_data['kx_um'], slice_data['ky_um'], 
                          slice_data['Rs'], shading='auto', cmap='viridis')
            plt.colorbar(label='Rs')
            plt.title(f'Rs at {energy:.1f} eV')
            plt.xlabel('kx (μm⁻¹)')
            plt.ylabel('ky (μm⁻¹)')
            
            print(f"✓ Extracted 2D slice at {energy:.1f} eV")
        
        plt.tight_layout()
        plt.savefig('jones_unified_2d_slices.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✓ 2D slice extraction completed")
        
    except Exception as e:
        print(f"✗ 2D slice extraction failed: {str(e)}")

def demo_1d_dispersion_extraction():
    """Demo 4: 1D dispersion extraction from 3D results"""
    print("\n" + "="*60)
    print("DEMO 4: 1D Dispersion Extraction")
    print("="*60)
    
    # Create structure
    layers = create_basic_structure()
    analyzer = Jones3DAnalyzer(layers)
    
    # Quick 3D calculation for dispersion extraction
    params = CalculationParams(
        energy_start_eV=1.0,
        energy_stop_eV=3.0,
        energy_points=31,
        kx_start_um=0.0,
        kx_stop_um=25.0,
        kx_points=26,
        ky_start_um=0.0,
        ky_stop_um=25.0,
        ky_points=26,
        es_amplitude=1.0,
        ep_amplitude=0.0,
        phase_diff_deg=0.0,
        resolution=GridResolution.MEDIUM
    )
    
    print("Running 3D calculation for dispersion extraction...")
    try:
        results = analyzer.calculate_3d_full(params)
        
        # Extract 1D dispersion along different k-directions
        kx_values = [0.0, 5.0, 10.0]
        
        plt.figure(figsize=(15, 5))
        
        for i, kx_val in enumerate(kx_values):
            disp_data = analyzer.get_1d_dispersion(results, 'kx', kx_val)
            
            plt.subplot(1, 3, i+1)
            plt.pcolormesh(disp_data['k_um'], disp_data['energy_eV'], 
                          disp_data['Rs'], shading='auto', cmap='viridis')
            plt.colorbar(label='Rs')
            plt.title(f'Rs Dispersion (kx={kx_val:.1f} μm⁻¹)')
            plt.xlabel('ky (μm⁻¹)')
            plt.ylabel('Energy (eV)')
            
            print(f"✓ Extracted 1D dispersion at kx={kx_val:.1f} μm⁻¹")
        
        plt.tight_layout()
        plt.savefig('jones_unified_1d_dispersion.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✓ 1D dispersion extraction completed")
        
    except Exception as e:
        print(f"✗ 1D dispersion extraction failed: {str(e)}")

def plot_preview_results(results, title):
    """Plot preview results (simple 2x2 subplot)"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot Rs and Rp at middle energy
    mid_energy = results.energy_eV.shape[0] // 2
    energy_val = results.energy_eV[mid_energy]
    
    # Rs
    im1 = axes[0, 0].pcolormesh(results.kx_um, results.ky_um, 
                               results.Rs[mid_energy, :, :], 
                               shading='auto', cmap='viridis')
    axes[0, 0].set_title(f'Rs at {energy_val:.2f} eV')
    axes[0, 0].set_xlabel('kx (μm⁻¹)')
    axes[0, 0].set_ylabel('ky (μm⁻¹)')
    fig.colorbar(im1, ax=axes[0, 0])
    
    # Rp
    im2 = axes[0, 1].pcolormesh(results.kx_um, results.ky_um, 
                               results.Rp[mid_energy, :, :], 
                               shading='auto', cmap='viridis')
    axes[0, 1].set_title(f'Rp at {energy_val:.2f} eV')
    axes[0, 1].set_xlabel('kx (μm⁻¹)')
    axes[0, 1].set_ylabel('ky (μm⁻¹)')
    fig.colorbar(im2, ax=axes[0, 1])
    
    # S0
    im3 = axes[1, 0].pcolormesh(results.kx_um, results.ky_um, 
                               results.S0[mid_energy, :, :], 
                               shading='auto', cmap='viridis')
    axes[1, 0].set_title(f'S0 at {energy_val:.2f} eV')
    axes[1, 0].set_xlabel('kx (μm⁻¹)')
    axes[1, 0].set_ylabel('ky (μm⁻¹)')
    fig.colorbar(im3, ax=axes[1, 0])
    
    # S1
    im4 = axes[1, 1].pcolormesh(results.kx_um, results.ky_um, 
                               results.S1[mid_energy, :, :], 
                               shading='auto', cmap='RdBu_r', vmin=-1, vmax=1)
    axes[1, 1].set_title(f'S1 at {energy_val:.2f} eV')
    axes[1, 1].set_xlabel('kx (μm⁻¹)')
    axes[1, 1].set_ylabel('ky (μm⁻¹)')
    fig.colorbar(im4, ax=axes[1, 1])
    
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'jones_unified_{title.lower().replace(" ", "_")}.png', 
                dpi=300, bbox_inches='tight')
    plt.show()

def plot_detailed_results(results, title):
    """Plot detailed results (full 2x4 subplot)"""
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    
    # Plot all 8 quantities at middle energy
    mid_energy = results.energy_eV.shape[0] // 2
    energy_val = results.energy_eV[mid_energy]
    
    quantities = ['Rs', 'Rp', 'S0', 'S1', 'S2', 'S3', 'phi', 'chi']
    data_arrays = [results.Rs, results.Rp, results.S0, results.S1, 
                   results.S2, results.S3, results.phi, results.chi]
    
    # Color map and range settings
    cmaps = ['viridis', 'viridis', 'viridis', 'RdBu_r', 'RdBu_r', 'RdBu_r', 'hsv', 'RdBu_r']
    vmin_vals = [0, 0, 0, -1, -1, -1, -np.pi/2, -np.pi/4]
    vmax_vals = [1, 1, 2, 1, 1, 1, np.pi/2, np.pi/4]
    
    for i, (ax, title_str, data, cmap, vmin, vmax) in enumerate(zip(
        axes.flat, quantities, data_arrays, cmaps, vmin_vals, vmax_vals)):
        
        im = ax.pcolormesh(results.kx_um, results.ky_um, 
                          data[mid_energy, :, :], 
                          shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(f'{title_str} at {energy_val:.2f} eV')
        ax.set_xlabel('kx (μm⁻¹)')
        ax.set_ylabel('ky (μm⁻¹)')
        fig.colorbar(im, ax=ax)
    
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'jones_unified_{title.lower().replace(" ", "_").replace("-", "_")}.png', 
                dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """Main demo function"""
    print("="*60)
    print("UNIFIED JONES CALCULATOR DEMO")
    print("="*60)
    print("This demo showcases the unified Jones calculator with:")
    print("• Adaptive grid resolution (coarse/fine)")
    print("• 3D E-kx-ky analysis")
    print("• Memory-efficient calculations")
    print("• Flexible data extraction")
    print("="*60)
    
    # Run demos
    try:
        # Demo 1: Coarse grid preview
        coarse_results = demo_coarse_grid_preview()
        
        # Demo 2: Fine grid analysis
        fine_results = demo_fine_grid_analysis()
        
        # Demo 3: 2D slice extraction
        demo_2d_slice_extraction()
        
        # Demo 4: 1D dispersion extraction
        demo_1d_dispersion_extraction()
        
        print("\n" + "="*60)
        print("✓ ALL DEMOS COMPLETED SUCCESSFULLY!")
        print("="*60)
        print("Generated files:")
        print("• jones_unified_coarse_grid_preview.png")
        print("• jones_unified_fine_grid_analysis_plasmonic_structure.png")
        print("• jones_unified_2d_slices.png")
        print("• jones_unified_1d_dispersion.png")
        print("="*60)
        
    except Exception as e:
        print(f"\n✗ Demo failed: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 