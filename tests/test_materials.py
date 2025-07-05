#!/usr/bin/env python3
"""
Tests for TMM materials module.

This module contains unit tests for the material system,
including Material class, MaterialLoader, and built-in materials.
"""

import pytest
import numpy as np
import tempfile
from pathlib import Path

from tmm.materials.material import Material
from tmm.materials.builtin_materials import (
    get_builtin_material, 
    list_builtin_materials,
    get_material_info
)
from tmm.materials.material_loader import MaterialLoader

class TestMaterial:
    """Test cases for the Material class."""
    
    def test_material_creation_single_energy(self):
        """Test creating a material with single energy point."""
        energy = 2.0
        eps = 2.25 + 0.1j
        
        material = Material("Test", energy, eps)
        
        assert material.name == "Test"
        assert material.energy_eV == energy
        assert material.eps_xx == eps
        assert not material.is_anisotropic
    
    def test_material_creation_array(self):
        """Test creating a material with energy array."""
        energy = np.linspace(1.0, 3.0, 10)
        eps = 2.0 + 0.1j * np.ones_like(energy)
        
        material = Material("Test", energy, eps)
        
        assert material.name == "Test"
        assert len(material.energy_eV) == 10
        assert not material.is_anisotropic
    
    def test_anisotropic_material(self):
        """Test creating an anisotropic material."""
        energy = 2.0
        eps_x = 2.0 + 0.1j
        eps_y = 2.5 + 0.1j
        eps_z = 3.0 + 0.1j
        
        material = Material("Anisotropic", energy, eps_x, eps_y, eps_z)
        
        assert material.is_anisotropic
        assert material.eps_xx == eps_x
        assert material.eps_yy == eps_y
        assert material.eps_zz == eps_z
    
    def test_dielectric_tensor_interpolation(self):
        """Test dielectric tensor interpolation."""
        energy_base = np.array([1.0, 2.0, 3.0])
        eps_base = np.array([1.0, 2.0, 3.0]) + 0.1j
        
        material = Material("Test", energy_base, eps_base)
        
        # Test interpolation
        energy_interp = np.array([1.5, 2.5])
        eps_tensor = material.get_dielectric_tensor(energy_interp)
        
        assert eps_tensor.shape == (2, 3)
        assert np.allclose(eps_tensor[:, 0].real, [1.5, 2.5], atol=1e-10)
    
    def test_refractive_index_calculation(self):
        """Test refractive index calculation."""
        energy = 2.0
        n = 1.5
        k = 0.1
        eps = (n + 1j * k)**2
        
        material = Material("Test", energy, eps)
        n_k = material.get_refractive_index(energy)
        
        assert np.isclose(n_k['n_x'], n, atol=1e-10)
        assert np.isclose(n_k['k_x'], k, atol=1e-10)
    
    def test_material_serialization(self):
        """Test material to/from dictionary conversion."""
        energy = np.linspace(1.0, 3.0, 5)
        eps = 2.0 + 0.1j * np.ones_like(energy)
        
        material1 = Material("Test", energy, eps, description="Test material")
        
        # Convert to dictionary
        data = material1.to_dict()
        
        # Create from dictionary
        material2 = Material.from_dict(data)
        
        assert material2.name == material1.name
        assert material2.description == material1.description
        assert np.allclose(material2.energy_eV, material1.energy_eV)
        assert np.allclose(material2.eps_xx, material1.eps_xx)

class TestBuiltinMaterials:
    """Test cases for built-in materials."""
    
    def test_list_builtin_materials(self):
        """Test listing available built-in materials."""
        materials = list_builtin_materials()
        
        assert isinstance(materials, list)
        assert len(materials) > 0
        assert "Air" in materials
        assert "SiO2" in materials
    
    def test_air_material(self):
        """Test Air material creation."""
        energy = np.linspace(1.0, 3.0, 10)
        air = get_builtin_material("Air", energy)
        
        assert air.name == "Air"
        eps_tensor = air.get_dielectric_tensor(energy)
        
        # Air should have eps ≈ 1
        assert np.allclose(eps_tensor[:, 0], 1.0 + 0j, atol=1e-10)
    
    def test_sio2_material(self):
        """Test SiO2 material creation."""
        energy = 2.0
        sio2 = get_builtin_material("SiO2", energy, n=1.46)
        
        assert sio2.name == "SiO2"
        
        # Check refractive index
        n_k = sio2.get_refractive_index(energy)
        assert np.isclose(n_k['n_x'], 1.46, atol=1e-10)
        assert np.isclose(n_k['k_x'], 0.0, atol=1e-10)
    
    def test_unknown_material(self):
        """Test handling of unknown material."""
        with pytest.raises(ValueError, match="Unknown material"):
            get_builtin_material("NonexistentMaterial", 2.0)
    
    def test_material_info(self):
        """Test material information retrieval."""
        info = get_material_info("Air")
        assert isinstance(info, str)
        assert "air" in info.lower()

class TestMaterialLoader:
    """Test cases for MaterialLoader."""
    
    def test_wavelength_to_energy_conversion(self):
        """Test wavelength to energy conversion."""
        wavelength_um = np.array([0.5, 1.0, 2.0])  # micrometers
        energy_eV = MaterialLoader._wavelength_to_energy(wavelength_um, "um")
        
        # Check that E * λ ≈ hc
        hc_eV_um = 1.239842  # hc in eV·μm
        calculated_hc = energy_eV * wavelength_um
        
        assert np.allclose(calculated_hc, hc_eV_um, rtol=1e-4)
    
    def test_energy_unit_conversion(self):
        """Test energy unit conversions."""
        # Test eV (should be unchanged)
        energy_eV = np.array([1.0, 2.0, 3.0])
        converted = MaterialLoader._energy_to_eV(energy_eV, "eV")
        assert np.allclose(converted, energy_eV)
        
        # Test cm⁻¹ to eV
        wavenumber_cm = np.array([8065.54429])  # Should be 1 eV
        converted = MaterialLoader._energy_to_eV(wavenumber_cm, "cm-1")
        assert np.isclose(converted[0], 1.0, rtol=1e-6)
    
    def test_nk_file_loading(self):
        """Test loading material from n,k file."""
        # Create temporary n,k file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("# wavelength(nm) n k\n")
            f.write("500 1.5 0.01\n")
            f.write("600 1.48 0.008\n")
            f.write("700 1.46 0.006\n")
            temp_file = f.name
        
        try:
            material = MaterialLoader.load_from_nk_file(
                temp_file, 
                material_name="TestMaterial",
                wavelength_unit="nm"
            )
            
            assert material.name == "TestMaterial"
            assert len(material.energy_eV) == 3
            
            # Check that energy is in ascending order
            assert np.all(material.energy_eV[1:] > material.energy_eV[:-1])
            
        finally:
            Path(temp_file).unlink()
    
    def test_dielectric_file_loading(self):
        """Test loading material from dielectric function file."""
        # Create temporary dielectric file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("# energy(eV) eps_real eps_imag\n")
            f.write("1.0 2.0 0.1\n")
            f.write("2.0 2.5 0.2\n")
            f.write("3.0 3.0 0.3\n")
            temp_file = f.name
        
        try:
            material = MaterialLoader.load_from_dielectric_file(
                temp_file,
                energy_unit="eV"
            )
            
            assert len(material.energy_eV) == 3
            
            # Check first data point
            eps_tensor = material.get_dielectric_tensor(1.0)
            expected_eps = 2.0 + 0.1j
            assert np.isclose(eps_tensor[0], expected_eps)
            
        finally:
            Path(temp_file).unlink()
    
    def test_file_format_detection(self):
        """Test automatic file format detection."""
        # Create temporary files with different formats
        
        # n,k file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("# wavelength n k\n")
            f.write("500 1.5 0.01\n")
            nk_file = f.name
        
        # Dielectric file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("# energy eps_real eps_imag\n")
            f.write("1.0 2.0 0.1\n")
            eps_file = f.name
        
        try:
            assert MaterialLoader.detect_file_format(nk_file) == "nk"
            assert MaterialLoader.detect_file_format(eps_file) == "dielectric"
            
        finally:
            Path(nk_file).unlink()
            Path(eps_file).unlink()

def test_conversion_utility():
    """Test n,k to dielectric conversion utility."""
    # Create temporary n,k file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("500 1.5 0.01\n")
        f.write("600 1.48 0.008\n")
        input_file = f.name
    
    output_file = input_file.replace('.txt', '_converted.txt')
    
    try:
        MaterialLoader.convert_nk_to_dielectric(
            input_file, 
            output_file,
            wavelength_unit="nm"
        )
        
        # Check that output file was created
        assert Path(output_file).exists()
        
        # Load and check the converted data
        data = np.loadtxt(output_file)
        assert data.shape[1] == 7  # energy + 3 complex dielectric components
        
    finally:
        Path(input_file).unlink()
        if Path(output_file).exists():
            Path(output_file).unlink()

if __name__ == "__main__":
    # Run tests if script is executed directly
    pytest.main([__file__, "-v"]) 