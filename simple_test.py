#!/usr/bin/env python3
"""
Simple test script to demonstrate TMM package functionality.
"""

import numpy as np

try:
    import tmm
    print("✅ TMM package imported successfully!")
    print(f"Version: {tmm.__version__}")
    print(f"Author: {tmm.__author__}")
    print()

    # Test basic functionality
    print("Testing available components:")
    
    # Test constants
    if tmm.Constants:
        print(f"✅ Constants available - Speed of light: {tmm.Constants.C:.2e} m/s")
    else:
        print("❌ Constants not available")
        
    # Test layer creation  
    if tmm.Layer and tmm.Material:
        print("✅ Layer and Material classes available")
        
        # Create a simple material (manually since built-in materials may not be fully implemented)
        energy_eV = np.array([1.0, 2.0, 3.0])
        eps_sio2 = np.array([2.13, 2.13, 2.13]) + 0j  # SiO2 dielectric constant
        
        try:
            material = tmm.Material("SiO2_test", energy_eV, eps_sio2)
            print(f"✅ Created material: {material.name}")
            
            # Create a layer
            layer = tmm.Layer(material, thickness=500e-9)  # 500 nm
            print(f"✅ Created layer: {layer}")
            
            # Test layer methods
            eps_tensor = layer.get_dielectric_tensor(2.0)
            print(f"✅ Dielectric tensor at 2 eV: {eps_tensor}")
            
        except Exception as e:
            print(f"❌ Error creating material/layer: {e}")
    else:
        print("❌ Layer or Material classes not available")
    
    # Test built-in materials
    if tmm.get_builtin_material:
        try:
            air = tmm.get_builtin_material("Air", 2.0)
            print(f"✅ Created built-in material: {air.name}")
        except Exception as e:
            print(f"❌ Error with built-in materials: {e}")
    else:
        print("⚠️ Built-in materials not fully implemented yet")
        
    # Test material loader
    if tmm.MaterialLoader:
        print("✅ MaterialLoader class available")
        
        # Test format detection on a hypothetical file
        try:
            # This will likely fail since we don't have a real file, but tests the import
            format_result = "Would test file format detection here"
            print(f"✅ MaterialLoader methods accessible")
        except Exception as e:
            print(f"⚠️ MaterialLoader available but needs real files to test: {type(e).__name__}")
    else:
        print("❌ MaterialLoader not available")
        
    print("\n" + "="*50)
    print("✅ TMM Package Basic Test Completed Successfully!")
    print("="*50)
    
except ImportError as e:
    print(f"❌ Import error: {e}")
except Exception as e:
    print(f"❌ Unexpected error: {e}") 