"""
Material loader for reading material data from files.

This module provides functionality to load material optical properties
from various file formats commonly used in optical simulations.
"""

import numpy as np
from typing import Union, Optional, Dict, List, Tuple
from pathlib import Path
import json
from .material import Material
from ..utils.validation_utils import ValidationUtils

class MaterialLoader:
    """
    Loader for material data from various file formats.
    
    Supports loading material data from:
    - Text files with wavelength/energy, n, k columns
    - Text files with dielectric function data
    - JSON files with material definitions
    """
    
    @staticmethod
    def load_from_nk_file(
        file_path: Union[str, Path],
        material_name: Optional[str] = None,
        wavelength_unit: str = "um",
        skip_header: int = 0,
        delimiter: Optional[str] = None
    ) -> Material:
        """
        Load material from n,k data file.
        
        Parameters
        ----------
        file_path : str or Path
            Path to the data file
        material_name : str, optional
            Name of the material (defaults to filename)
        wavelength_unit : str, default="um"
            Unit of wavelength column ("um", "nm", "m")
        skip_header : int, default=0
            Number of header lines to skip
        delimiter : str, optional
            Column delimiter (auto-detected if None)
            
        Returns
        -------
        Material
            Material object with loaded data
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        # Load data
        data = np.loadtxt(file_path, skiprows=skip_header, delimiter=delimiter)
        
        if data.shape[1] < 3:
            raise ValueError("File must have at least 3 columns: wavelength, n, k")
        
        wavelength = data[:, 0]
        n = data[:, 1]
        k = data[:, 2]
        
        # Convert wavelength to energy
        energy_eV = MaterialLoader._wavelength_to_energy(wavelength, wavelength_unit)
        
        # Calculate dielectric function
        eps_complex = (n + 1j * k)**2
        
        # Sort by energy (ascending)
        sort_idx = np.argsort(energy_eV)
        energy_eV = energy_eV[sort_idx]
        eps_complex = eps_complex[sort_idx]
        
        # Create material
        if material_name is None:
            material_name = file_path.stem
        
        return Material(
            name=material_name,
            energy_eV=energy_eV,
            eps_xx=eps_complex,
            description=f"Material loaded from {file_path.name}",
            references=f"Data file: {file_path.name}"
        )
    
    @staticmethod
    def load_from_dielectric_file(
        file_path: Union[str, Path],
        material_name: Optional[str] = None,
        energy_unit: str = "eV",
        skip_header: int = 0,
        delimiter: Optional[str] = None
    ) -> Material:
        """
        Load material from dielectric function file.
        
        Expected format: energy/wavenumber, eps_real, eps_imag, [eps_y_real, eps_y_imag, eps_z_real, eps_z_imag]
        
        Parameters
        ----------
        file_path : str or Path
            Path to the data file
        material_name : str, optional
            Name of the material (defaults to filename)
        energy_unit : str, default="eV"
            Unit of energy column ("eV", "cm-1", "THz")
        skip_header : int, default=0
            Number of header lines to skip
        delimiter : str, optional
            Column delimiter (auto-detected if None)
            
        Returns
        -------
        Material
            Material object with loaded data
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        # Load data
        data = np.loadtxt(file_path, skiprows=skip_header, delimiter=delimiter)
        
        if data.shape[1] < 3:
            raise ValueError("File must have at least 3 columns: energy, eps_real, eps_imag")
        
        energy = data[:, 0]
        eps_x_real = data[:, 1]
        eps_x_imag = data[:, 2]
        
        # Convert energy to eV if needed
        energy_eV = MaterialLoader._energy_to_eV(energy, energy_unit)
        
        # Build dielectric function
        eps_x = eps_x_real + 1j * eps_x_imag
        
        # Check for anisotropic data
        if data.shape[1] >= 7:
            eps_y_real = data[:, 3]
            eps_y_imag = data[:, 4]
            eps_z_real = data[:, 5]
            eps_z_imag = data[:, 6]
            
            eps_y = eps_y_real + 1j * eps_y_imag
            eps_z = eps_z_real + 1j * eps_z_imag
        else:
            eps_y = eps_x.copy()
            eps_z = eps_x.copy()
        
        # Sort by energy (ascending)
        sort_idx = np.argsort(energy_eV)
        energy_eV = energy_eV[sort_idx]
        eps_x = eps_x[sort_idx]
        eps_y = eps_y[sort_idx]
        eps_z = eps_z[sort_idx]
        
        # Create material
        if material_name is None:
            material_name = file_path.stem
        
        return Material(
            name=material_name,
            energy_eV=energy_eV,
            eps_xx=eps_x,
            eps_yy=eps_y,
            eps_zz=eps_z,
            description=f"Material loaded from {file_path.name}",
            references=f"Data file: {file_path.name}"
        )
    
    @staticmethod
    def load_from_json(file_path: Union[str, Path]) -> Material:
        """
        Load material from JSON file.
        
        Parameters
        ----------
        file_path : str or Path
            Path to the JSON file
            
        Returns
        -------
        Material
            Material object with loaded data
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        return Material.from_dict(data)
    
    @staticmethod
    def save_to_json(material: Material, file_path: Union[str, Path]) -> None:
        """
        Save material to JSON file.
        
        Parameters
        ----------
        material : Material
            Material object to save
        file_path : str or Path
            Path to save the JSON file
        """
        file_path = Path(file_path)
        
        data = material.to_dict()
        
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    @staticmethod
    def detect_file_format(file_path: Union[str, Path]) -> str:
        """
        Detect the format of a material data file.
        
        Parameters
        ----------
        file_path : str or Path
            Path to the file
            
        Returns
        -------
        str
            Detected format: "nk", "dielectric", "json", "unknown"
        """
        file_path = Path(file_path)
        
        if file_path.suffix.lower() == '.json':
            return "json"
        
        if not file_path.exists():
            return "unknown"
        
        try:
            # Try to read first few lines
            with open(file_path, 'r') as f:
                lines = [f.readline().strip() for _ in range(5)]
            
            # Look for header information
            header_text = ' '.join(lines).lower()
            
            if 'wavelength' in header_text and ('n' in header_text or 'k' in header_text):
                return "nk"
            elif 'dielectric' in header_text or 'eps' in header_text:
                return "dielectric"
            elif 'energy' in header_text or 'wavenumber' in header_text:
                return "dielectric"
            else:
                # Try to infer from data structure
                data = np.loadtxt(file_path, max_rows=10)
                if data.shape[1] == 3:
                    # Could be wavelength, n, k
                    return "nk"
                elif data.shape[1] >= 3:
                    # Could be energy, eps_real, eps_imag
                    return "dielectric"
                else:
                    return "unknown"
        except Exception:
            return "unknown"
    
    @staticmethod
    def _wavelength_to_energy(wavelength: np.ndarray, unit: str) -> np.ndarray:
        """Convert wavelength to energy in eV."""
        # Convert to meters
        if unit == "um":
            wavelength_m = wavelength * 1e-6
        elif unit == "nm":
            wavelength_m = wavelength * 1e-9
        elif unit == "m":
            wavelength_m = wavelength
        else:
            raise ValueError(f"Unknown wavelength unit: {unit}")
        
        # Convert to energy in eV
        from ..utils.constants import Constants
        energy_J = Constants.H * Constants.C / wavelength_m
        energy_eV = energy_J / Constants.EV_TO_J
        
        return energy_eV
    
    @staticmethod
    def _energy_to_eV(energy: np.ndarray, unit: str) -> np.ndarray:
        """Convert energy to eV."""
        if unit == "eV":
            return energy
        elif unit == "cm-1":
            from ..utils.constants import Constants
            return energy * Constants.WAVENUMBER_TO_EV
        elif unit == "THz":
            from ..utils.constants import Constants
            return energy * Constants.H / Constants.EV_TO_J * 1e12
        else:
            raise ValueError(f"Unknown energy unit: {unit}")
    
    @staticmethod
    def convert_nk_to_dielectric(
        input_path: Union[str, Path],
        output_path: Union[str, Path],
        wavelength_unit: str = "um",
        skip_header: int = 0
    ) -> None:
        """
        Convert n,k file to dielectric function file.
        
        Parameters
        ----------
        input_path : str or Path
            Path to n,k input file
        output_path : str or Path
            Path to dielectric output file
        wavelength_unit : str, default="um"
            Unit of wavelength in input file
        skip_header : int, default=0
            Number of header lines to skip
        """
        input_path = Path(input_path)
        output_path = Path(output_path)
        
        # Load n,k data
        data = np.loadtxt(input_path, skiprows=skip_header)
        wavelength = data[:, 0]
        n = data[:, 1]
        k = data[:, 2]
        
        # Convert to energy and dielectric function
        energy_eV = MaterialLoader._wavelength_to_energy(wavelength, wavelength_unit)
        eps = (n + 1j * k)**2
        
        # Sort by energy
        sort_idx = np.argsort(energy_eV)
        energy_eV = energy_eV[sort_idx]
        eps = eps[sort_idx]
        
        # Save dielectric function data
        header = "Energy(eV) eps_real eps_imag eps_y_real eps_y_imag eps_z_real eps_z_imag"
        output_data = np.column_stack([
            energy_eV,
            eps.real, eps.imag,
            eps.real, eps.imag,  # duplicate for y
            eps.real, eps.imag   # duplicate for z
        ])
        
        np.savetxt(output_path, output_data, header=header, fmt="%.6f")
        
        print(f"Converted {input_path.name} -> {output_path.name}")

def load_material_database(directory: Union[str, Path]) -> Dict[str, Material]:
    """
    Load all materials from a directory.
    
    Parameters
    ----------
    directory : str or Path
        Directory containing material files
        
    Returns
    -------
    Dict[str, Material]
        Dictionary mapping material names to Material objects
    """
    directory = Path(directory)
    
    if not directory.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    materials = {}
    
    for file_path in directory.glob("*"):
        if file_path.is_file():
            try:
                file_format = MaterialLoader.detect_file_format(file_path)
                
                if file_format == "nk":
                    material = MaterialLoader.load_from_nk_file(file_path)
                elif file_format == "dielectric":
                    material = MaterialLoader.load_from_dielectric_file(file_path)
                elif file_format == "json":
                    material = MaterialLoader.load_from_json(file_path)
                else:
                    continue
                
                materials[material.name] = material
                print(f"Loaded material: {material.name}")
                
            except Exception as e:
                print(f"Warning: Could not load {file_path.name}: {e}")
    
    return materials 