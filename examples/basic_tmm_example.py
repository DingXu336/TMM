#!/usr/bin/env python3
"""
Basic TMM Example: Simple Fabry-Pérot Cavity

This example demonstrates basic usage of the TMM package for calculating
the reflectivity of a simple Fabry-Pérot cavity structure.

Author: Ding Xu
"""

import numpy as np
import matplotlib.pyplot as plt

# Import TMM package components
from tmm.materials import get_builtin_material
from tmm.core import Layer, create_substrate, create_superstrate, create_film

def create_fabry_perot_structure():
    """
    Create a simple Fabry-Pérot cavity structure:
    Air / SiO2 (500 nm) / Air
    """
    print("Creating Fabry-Pérot cavity structure...")
    
    # Define energy range for calculations
    energy_eV = np.linspace(1.5, 3.0, 200)
    
    # Create materials
    air = get_builtin_material("Air", energy_eV)
    sio2 = get_builtin_material("SiO2", energy_eV, n=1.46)
    
    # Create layer structure
    layers = [
        create_superstrate(air, "Air_top"),           # Semi-infinite air (top)
        create_film(sio2, thickness_nm=500, "SiO2_cavity"),  # 500 nm SiO2 cavity
        create_substrate(air, "Air_bottom")           # Semi-infinite air (bottom)
    ]
    
    print(f"Created {len(layers)} layers:")
    for i, layer in enumerate(layers):
        print(f"  {i+1}. {layer}")
    
    return layers, energy_eV

def calculate_basic_reflectivity(layers, energy_eV):
    """
    Calculate basic reflectivity spectrum.
    
    Note: This is a simplified example. The full TMMCalculator 
    implementation would be more comprehensive.
    """
    print("\nCalculating reflectivity spectrum...")
    
    # For this example, we'll create a simplified calculation
    # In the full implementation, this would use the TMMCalculator class
    
    # Get material properties
    air_eps = layers[0].get_dielectric_tensor(energy_eV)
    sio2_eps = layers[1].get_dielectric_tensor(energy_eV)
    
    # Simple Fabry-Pérot approximation for demonstration
    # In reality, the TMM calculation would be more rigorous
    
    wavelength_vacuum = np.array([6.626e-34 * 3e8 / (E * 1.602e-19) for E in energy_eV])
    n_sio2 = np.sqrt(sio2_eps[:, 0].real)  # Refractive index
    
    # Optical thickness
    thickness = layers[1].thickness
    optical_thickness = 2 * n_sio2 * thickness / wavelength_vacuum
    
    # Simple Fabry-Pérot reflectivity (approximate)
    # R = |r|² where r is the reflection coefficient
    # This is simplified - full TMM would handle multiple reflections properly
    r12 = (1 - n_sio2) / (1 + n_sio2)  # Air-SiO2 interface
    r21 = -r12  # SiO2-Air interface
    
    # Include multiple reflections (simplified)
    phase = 2 * np.pi * optical_thickness
    numerator = r12 + r21 * np.exp(2j * phase)
    denominator = 1 + r12 * r21 * np.exp(2j * phase)
    
    r_total = numerator / denominator
    reflectivity = np.abs(r_total)**2
    
    return reflectivity

def plot_results(energy_eV, reflectivity):
    """Plot the reflectivity spectrum."""
    print("\nPlotting results...")
    
    # Convert energy to wavelength for additional x-axis
    wavelength_nm = 1239.84 / energy_eV  # hc/E in nm
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot vs energy
    ax1.plot(energy_eV, reflectivity, 'b-', linewidth=2, label='Reflectivity')
    ax1.set_xlabel('Energy (eV)')
    ax1.set_ylabel('Reflectivity')
    ax1.set_title('Fabry-Pérot Cavity Reflectivity Spectrum')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim(0, 1)
    
    # Plot vs wavelength
    ax2.plot(wavelength_nm, reflectivity, 'r-', linewidth=2, label='Reflectivity')
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel('Reflectivity')
    ax2.set_title('Fabry-Pérot Cavity Reflectivity Spectrum')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_ylim(0, 1)
    ax2.invert_xaxis()  # Higher energy = shorter wavelength
    
    plt.tight_layout()
    plt.savefig('fabry_perot_reflectivity.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Plot saved as 'fabry_perot_reflectivity.png'")

def analyze_cavity_modes(energy_eV, reflectivity):
    """Analyze cavity modes from the reflectivity spectrum."""
    print("\nAnalyzing cavity modes...")
    
    # Find minima in reflectivity (transmission maxima)
    from scipy.signal import find_peaks
    
    # Invert reflectivity to find transmission peaks
    transmission = 1 - reflectivity
    peaks, properties = find_peaks(transmission, height=0.5, distance=10)
    
    if len(peaks) > 0:
        print(f"Found {len(peaks)} cavity modes:")
        for i, peak_idx in enumerate(peaks):
            energy_mode = energy_eV[peak_idx]
            wavelength_mode = 1239.84 / energy_mode
            transmission_max = transmission[peak_idx]
            print(f"  Mode {i+1}: E = {energy_mode:.3f} eV, λ = {wavelength_mode:.1f} nm, T = {transmission_max:.3f}")
    else:
        print("No clear cavity modes found in this energy range")

def main():
    """Main example function."""
    print("=" * 60)
    print("TMM Package - Basic Fabry-Pérot Cavity Example")
    print("=" * 60)
    
    try:
        # Create structure
        layers, energy_eV = create_fabry_perot_structure()
        
        # Calculate reflectivity
        reflectivity = calculate_basic_reflectivity(layers, energy_eV)
        
        # Analyze results
        analyze_cavity_modes(energy_eV, reflectivity)
        
        # Plot results
        plot_results(energy_eV, reflectivity)
        
        print("\n" + "=" * 60)
        print("Example completed successfully!")
        print("=" * 60)
        
    except ImportError as e:
        print(f"Import error: {e}")
        print("Please install required dependencies:")
        print("  pip install numpy matplotlib scipy")
    except Exception as e:
        print(f"Error: {e}")
        print("Please check your installation and try again.")

if __name__ == "__main__":
    main() 