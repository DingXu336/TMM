#!/usr/bin/env python3
"""
Working TMM Example - Demonstrates Implemented Functionality

This example shows how to use the TMM package with the components that are 
currently implemented and working.

Author: Ding Xu
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    print("=" * 60)
    print("TMM Package - Working Example")
    print("=" * 60)
    
    try:
        # Import TMM package
        import tmm
        print(f"✅ TMM package loaded successfully!")
        print(f"   Version: {tmm.__version__}")
        print(f"   Author: {tmm.__author__}")
        print()
        
        # 1. Create materials manually
        print("1. Creating Materials")
        print("-" * 20)
        
        # Energy range for calculations
        energy_eV = np.linspace(1.5, 3.5, 50)
        print(f"   Energy range: {energy_eV[0]:.1f} - {energy_eV[-1]:.1f} eV ({len(energy_eV)} points)")
        
        # Create Air (vacuum)
        eps_air = np.ones_like(energy_eV, dtype=complex)  # n ≈ 1 for air
        air = tmm.Material("Air", energy_eV, eps_air)
        print(f"   ✅ Created: {air}")
        
        # Create SiO2 (glass)
        n_sio2 = 1.46  # Typical refractive index for SiO2
        eps_sio2 = (n_sio2**2) * np.ones_like(energy_eV, dtype=complex)  
        sio2 = tmm.Material("SiO2", energy_eV, eps_sio2)
        print(f"   ✅ Created: {sio2}")
        
        # Create Gold (with simple Drude model approximation)
        # This is a simplified model - real gold would have more complex dispersion
        eps_gold_real = np.ones_like(energy_eV) * (-10)  # Simplified: mostly negative real part
        eps_gold_imag = np.ones_like(energy_eV) * 2     # Some absorption
        eps_gold = eps_gold_real + 1j * eps_gold_imag
        gold = tmm.Material("Gold_simplified", energy_eV, eps_gold)
        print(f"   ✅ Created: {gold}")
        print()
        
        # 2. Create layer structure
        print("2. Creating Layer Structure")
        print("-" * 20)
        
        # Create a simple Fabry-Pérot cavity: Air / SiO2 (500nm) / Air
        layers = [
            tmm.create_superstrate(air, "Air_top"),          # Semi-infinite air above
                         tmm.create_film(sio2, 500, "SiO2_cavity"),  # 500 nm SiO2 cavity
            tmm.create_substrate(air, "Air_bottom")          # Semi-infinite air below
        ]
        
        print(f"   Layer structure ({len(layers)} layers):")
        for i, layer in enumerate(layers):
            print(f"     {i+1}. {layer}")
        print()
        
        # 3. Analyze material properties
        print("3. Material Properties Analysis")
        print("-" * 20)
        
        # Test energy point
        test_energy = 2.0  # eV
        
        for material in [air, sio2, gold]:
            # Get optical constants
            optical_data = material.get_optical_constants(test_energy)
            n = optical_data['n_x']
            k = optical_data['k_x']
            eps = optical_data['eps_xx']
            
            print(f"   {material.name} at {test_energy} eV:")
            print(f"     n = {n:.3f}, k = {k:.3f}")
            print(f"     ε = {eps:.3f}")
        print()
        
        # 4. Calculate some basic optical properties
        print("4. Basic Optical Calculations")
        print("-" * 20)
        
        # Calculate optical thickness of the SiO2 layer
        sio2_layer = layers[1]  # The middle layer
        optical_thickness = sio2_layer.get_optical_thickness(test_energy)
        
        print(f"   SiO2 layer optical thickness at {test_energy} eV:")
        print(f"     Physical thickness: {sio2_layer.thickness_nm:.0f} nm")
        print(f"     Optical thickness: {optical_thickness:.3f} wavelengths")
        
        # Check if materials are absorbing
        for material in [air, sio2, gold]:
            is_absorbing = material.get_optical_constants(test_energy)['k_x'] > 1e-6
            status = "absorbing" if is_absorbing else "transparent"
            print(f"     {material.name}: {status}")
        print()
        
        # 5. Plot material dispersion
        print("5. Plotting Material Dispersion")
        print("-" * 20)
        
        # Calculate optical constants over energy range
        materials_data = {}
        for material in [air, sio2, gold]:
            data = material.get_optical_constants(energy_eV)
            materials_data[material.name] = {
                'n': data['n_x'],
                'k': data['k_x'],
                'eps_real': data['eps_xx'].real,
                'eps_imag': data['eps_xx'].imag
            }
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle('Material Optical Properties vs Energy', fontsize=14)
        
        colors = {'Air': 'blue', 'SiO2': 'green', 'Gold_simplified': 'gold'}
        
        for name, data in materials_data.items():
            color = colors.get(name, 'black')
            
            # Refractive index
            axes[0,0].plot(energy_eV, data['n'], color=color, label=name, linewidth=2)
            axes[0,0].set_ylabel('Refractive index (n)')
            axes[0,0].set_title('Refractive Index')
            axes[0,0].legend()
            axes[0,0].grid(True, alpha=0.3)
            
            # Extinction coefficient  
            axes[0,1].semilogy(energy_eV, np.maximum(data['k'], 1e-6), 
                              color=color, label=name, linewidth=2)
            axes[0,1].set_ylabel('Extinction coefficient (k)')
            axes[0,1].set_title('Extinction Coefficient (log scale)')
            axes[0,1].legend()
            axes[0,1].grid(True, alpha=0.3)
            
            # Dielectric function real part
            axes[1,0].plot(energy_eV, data['eps_real'], color=color, label=name, linewidth=2)
            axes[1,0].set_xlabel('Energy (eV)')
            axes[1,0].set_ylabel('Re(ε)')
            axes[1,0].set_title('Dielectric Function (Real)')
            axes[1,0].legend()
            axes[1,0].grid(True, alpha=0.3)
            
            # Dielectric function imaginary part
            axes[1,1].plot(energy_eV, data['eps_imag'], color=color, label=name, linewidth=2)
            axes[1,1].set_xlabel('Energy (eV)')
            axes[1,1].set_ylabel('Im(ε)')
            axes[1,1].set_title('Dielectric Function (Imaginary)')
            axes[1,1].legend()
            axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('tmm_materials_demo.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"   ✅ Plot saved as 'tmm_materials_demo.png'")
        print()
        
        # 6. Summary
        print("6. Summary")
        print("-" * 20)
        print("   ✅ Successfully demonstrated:")
        print("     • Material creation and management")
        print("     • Layer structure definition")  
        print("     • Optical constants calculation")
        print("     • Basic optical properties")
        print("     • Material dispersion analysis")
        print("     • Data visualization")
        print()
        print("   📋 Next steps for full implementation:")
        print("     • Complete TMM calculation engine")
        print("     • Reflectivity/transmission calculations")
        print("     • Jones matrix analysis") 
        print("     • GUI applications")
        
        print("\n" + "=" * 60)
        print("✅ TMM Working Example Completed Successfully!")
        print("=" * 60)
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("Please ensure the TMM package is installed: python3 -m pip install -e .")
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 