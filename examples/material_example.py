#!/usr/bin/env python3
"""
Material System Example

This example demonstrates how to work with materials in the TMM package:
- Using built-in materials
- Loading custom materials from files
- Creating and manipulating material properties

Author: Ding Xu
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Import TMM material components
from tmm.materials import (
    get_builtin_material, 
    list_builtin_materials,
    get_material_info,
    Material,
    MaterialLoader
)

def demonstrate_builtin_materials():
    """Demonstrate usage of built-in materials."""
    print("Built-in Materials Demo")
    print("-" * 40)
    
    # List all available built-in materials
    available_materials = list_builtin_materials()
    print(f"Available built-in materials: {', '.join(available_materials)}")
    
    # Energy range for calculations
    energy_eV = np.linspace(0.5, 4.0, 100)
    
    # Create different materials
    materials = {}
    
    # Air
    materials['Air'] = get_builtin_material("Air", energy_eV)
    print(f"\nCreated: {materials['Air']}")
    
    # SiO2 with default refractive index
    materials['SiO2'] = get_builtin_material("SiO2", energy_eV)
    print(f"Created: {materials['SiO2']}")
    
    # SiO2 with custom refractive index
    materials['SiO2_custom'] = get_builtin_material("SiO2", energy_eV, n=1.5)
    print(f"Created: {materials['SiO2_custom']}")
    
    # Silicon with dispersion
    materials['Si'] = get_builtin_material("Si", energy_eV, use_dispersion=True)
    print(f"Created: {materials['Si']}")
    
    # Gold with Drude model
    materials['Au'] = get_builtin_material("Au", energy_eV)
    print(f"Created: {materials['Au']}")
    
    return materials, energy_eV

def plot_material_properties(materials, energy_eV):
    """Plot optical properties of different materials."""
    print("\nPlotting material properties...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Optical Properties of Different Materials', fontsize=16)
    
    # Colors for different materials
    colors = {'Air': 'blue', 'SiO2': 'green', 'SiO2_custom': 'lightgreen', 
              'Si': 'red', 'Au': 'gold'}
    
    for name, material in materials.items():
        if name == 'SiO2_custom':  # Skip to avoid clutter
            continue
            
        color = colors.get(name, 'black')
        
        # Get optical constants
        optical_data = material.get_optical_constants(energy_eV)
        
        # Dielectric function real part
        ax = axes[0, 0]
        ax.plot(energy_eV, optical_data['eps_xx'].real, 
                color=color, label=f'{name} Re(ε)', linewidth=2)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Re(ε)')
        ax.set_title('Dielectric Function (Real Part)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Dielectric function imaginary part
        ax = axes[0, 1]
        ax.plot(energy_eV, optical_data['eps_xx'].imag, 
                color=color, label=f'{name} Im(ε)', linewidth=2)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Im(ε)')
        ax.set_title('Dielectric Function (Imaginary Part)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Refractive index
        ax = axes[1, 0]
        ax.plot(energy_eV, optical_data['n_x'], 
                color=color, label=f'{name} n', linewidth=2)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Refractive Index (n)')
        ax.set_title('Refractive Index')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Extinction coefficient
        ax = axes[1, 1]
        # Use log scale for extinction coefficient
        ax.semilogy(energy_eV, np.maximum(optical_data['k_x'], 1e-6), 
                   color=color, label=f'{name} k', linewidth=2)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Extinction Coefficient (k)')
        ax.set_title('Extinction Coefficient (log scale)')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('material_properties.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Material properties plot saved as 'material_properties.png'")

def create_custom_material():
    """Create a custom material from data."""
    print("\nCreating Custom Material")
    print("-" * 40)
    
    # Energy range
    energy_eV = np.linspace(1.0, 3.0, 50)
    
    # Create synthetic material data (example: dispersive material)
    # Lorentzian oscillator model
    E0 = 2.0  # Resonance energy (eV)
    gamma = 0.1  # Damping (eV)
    f = 1.0  # Oscillator strength
    eps_inf = 2.0  # High-frequency dielectric constant
    
    # Calculate dielectric function
    omega = energy_eV
    omega0 = E0
    eps_custom = eps_inf + f * omega0**2 / (omega0**2 - omega**2 - 1j * gamma * omega)
    
    # Create material
    custom_material = Material(
        name="Custom_Lorentzian",
        energy_eV=energy_eV,
        eps_xx=eps_custom,
        description="Custom material with Lorentzian dispersion",
        references="Created for TMM example"
    )
    
    print(f"Created custom material: {custom_material}")
    
    # Get and display optical constants
    optical_data = custom_material.get_optical_constants(energy_eV)
    
    print(f"Energy range: {energy_eV[0]:.2f} - {energy_eV[-1]:.2f} eV")
    print(f"Refractive index range: {optical_data['n_x'].min():.3f} - {optical_data['n_x'].max():.3f}")
    print(f"Extinction coefficient range: {optical_data['k_x'].min():.3f} - {optical_data['k_x'].max():.3f}")
    
    return custom_material

def demonstrate_material_loading():
    """Demonstrate material loading from files."""
    print("\nMaterial Loading Demo")
    print("-" * 40)
    
    # Create example data files
    create_example_material_files()
    
    try:
        # Load from n,k file
        print("Loading material from n,k file...")
        material_nk = MaterialLoader.load_from_nk_file(
            "example_material_nk.txt",
            material_name="Example_nk",
            wavelength_unit="nm"
        )
        print(f"Loaded: {material_nk}")
        
        # Load from dielectric file
        print("Loading material from dielectric file...")
        material_eps = MaterialLoader.load_from_dielectric_file(
            "example_material_eps.txt",
            material_name="Example_eps",
            energy_unit="eV"
        )
        print(f"Loaded: {material_eps}")
        
        # Save material to JSON
        print("Saving material to JSON...")
        MaterialLoader.save_to_json(material_nk, "example_material.json")
        print("Material saved to 'example_material.json'")
        
        # Load from JSON
        print("Loading material from JSON...")
        material_json = MaterialLoader.load_from_json("example_material.json")
        print(f"Loaded from JSON: {material_json}")
        
        return [material_nk, material_eps, material_json]
        
    except Exception as e:
        print(f"Error loading materials: {e}")
        return []

def create_example_material_files():
    """Create example material data files for demonstration."""
    
    # Create n,k data file
    wavelength_nm = np.linspace(400, 800, 20)
    n_data = 1.5 + 0.01 * np.sin(wavelength_nm / 100)  # Slight dispersion
    k_data = 0.001 * np.ones_like(wavelength_nm)  # Small absorption
    
    nk_data = np.column_stack([wavelength_nm, n_data, k_data])
    np.savetxt("example_material_nk.txt", nk_data, 
               header="wavelength(nm) n k", fmt="%.6f")
    
    # Create dielectric function file
    energy_eV = np.linspace(1.5, 3.1, 20)
    eps_real = 2.25 + 0.1 * np.sin(energy_eV)  # Slight dispersion
    eps_imag = 0.01 * np.ones_like(energy_eV)  # Small absorption
    
    eps_data = np.column_stack([energy_eV, eps_real, eps_imag, 
                               eps_real, eps_imag, eps_real, eps_imag])
    np.savetxt("example_material_eps.txt", eps_data,
               header="energy(eV) eps_x_real eps_x_imag eps_y_real eps_y_imag eps_z_real eps_z_imag",
               fmt="%.6f")

def main():
    """Main example function."""
    print("=" * 60)
    print("TMM Package - Material System Example")
    print("=" * 60)
    
    try:
        # Demonstrate built-in materials
        materials, energy_eV = demonstrate_builtin_materials()
        
        # Plot material properties
        plot_material_properties(materials, energy_eV)
        
        # Create custom material
        custom_material = create_custom_material()
        
        # Add custom material to collection
        materials['Custom'] = custom_material
        
        # Demonstrate material loading
        loaded_materials = demonstrate_material_loading()
        
        # Display material information
        print("\nMaterial Information:")
        print("-" * 40)
        for name in ['Air', 'SiO2', 'Si', 'Au']:
            info = get_material_info(name)
            print(f"{name}: {info}")
        
        print("\n" + "=" * 60)
        print("Material example completed successfully!")
        print("=" * 60)
        
    except ImportError as e:
        print(f"Import error: {e}")
        print("Please install required dependencies:")
        print("  pip install numpy matplotlib")
    except Exception as e:
        print(f"Error: {e}")
        print("Please check your installation and try again.")

if __name__ == "__main__":
    main() 