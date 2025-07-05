"""
Built-in materials for TMM calculations.

This module provides commonly used materials with their optical properties.
These materials can be used directly in TMM calculations without loading
external data files.
"""

import numpy as np
from typing import Union, Optional
from .material import Material

def get_air(energy_eV: Union[float, np.ndarray]) -> Material:
    """
    Create Air material with refractive index ≈ 1.
    
    Parameters
    ----------
    energy_eV : float or np.ndarray
        Energy values in eV
        
    Returns
    -------
    Material
        Air material object
    """
    energy_array = np.asarray(energy_eV)
    eps_air = np.ones_like(energy_array, dtype=complex)
    
    return Material(
        name="Air",
        energy_eV=energy_array,
        eps_xx=eps_air,
        description="Air with refractive index n ≈ 1.0",
        references="Standard air at room temperature and pressure"
    )

def get_vacuum(energy_eV: Union[float, np.ndarray]) -> Material:
    """
    Create Vacuum material with refractive index = 1.
    
    Parameters
    ----------
    energy_eV : float or np.ndarray
        Energy values in eV
        
    Returns
    -------
    Material
        Vacuum material object
    """
    energy_array = np.asarray(energy_eV)
    eps_vacuum = np.ones_like(energy_array, dtype=complex)
    
    return Material(
        name="Vacuum",
        energy_eV=energy_array,
        eps_xx=eps_vacuum,
        description="Perfect vacuum with refractive index n = 1.0",
        references="Theoretical perfect vacuum"
    )

def get_sio2(energy_eV: Union[float, np.ndarray], n: float = 1.46) -> Material:
    """
    Create SiO₂ (silica) material with constant refractive index.
    
    Parameters
    ----------
    energy_eV : float or np.ndarray
        Energy values in eV
    n : float, default=1.46
        Refractive index of SiO₂
        
    Returns
    -------
    Material
        SiO₂ material object
    """
    energy_array = np.asarray(energy_eV)
    eps_sio2 = np.full_like(energy_array, n**2, dtype=complex)
    
    return Material(
        name="SiO2",
        energy_eV=energy_array,
        eps_xx=eps_sio2,
        description=f"Silicon dioxide (SiO₂) with n = {n}",
        references="Typical value for fused silica glass"
    )

def get_si(energy_eV: Union[float, np.ndarray], 
           use_dispersion: bool = True) -> Material:
    """
    Create Silicon material with optional dispersion.
    
    Parameters
    ----------
    energy_eV : float or np.ndarray
        Energy values in eV
    use_dispersion : bool, default=True
        Whether to include dispersion (simplified model)
        
    Returns
    -------
    Material
        Silicon material object
    """
    energy_array = np.asarray(energy_eV)
    
    if use_dispersion:
        # Simplified dispersion model for Si
        # This is a very rough approximation - real Si has complex dispersion
        n_si = 3.5 + 0.1 / (energy_array**2 + 0.1)
        k_si = np.where(energy_array > 1.1, 0.01 * (energy_array - 1.1), 0)
        eps_si = (n_si + 1j * k_si)**2
    else:
        # Constant refractive index
        eps_si = np.full_like(energy_array, (3.5 + 0j)**2, dtype=complex)
    
    return Material(
        name="Si",
        energy_eV=energy_array,
        eps_xx=eps_si,
        description="Silicon with simplified dispersion model",
        references="Simplified model - use experimental data for accuracy"
    )

def get_au(energy_eV: Union[float, np.ndarray]) -> Material:
    """
    Create Gold material with Drude model.
    
    Parameters
    ----------
    energy_eV : float or np.ndarray
        Energy values in eV
        
    Returns
    -------
    Material
        Gold material object
    """
    energy_array = np.asarray(energy_eV)
    
    # Simplified Drude model for Au
    # Real Au has interband transitions - this is just a rough approximation
    omega_p = 8.9  # Plasma frequency in eV
    gamma = 0.07   # Damping parameter in eV
    
    eps_inf = 1.0  # High-frequency dielectric constant
    eps_au = eps_inf - (omega_p**2) / (energy_array**2 + 1j * gamma * energy_array)
    
    return Material(
        name="Au",
        energy_eV=energy_array,
        eps_xx=eps_au,
        description="Gold with simplified Drude model",
        references="Simplified Drude model - use experimental data for accuracy"
    )

def get_ag(energy_eV: Union[float, np.ndarray]) -> Material:
    """
    Create Silver material with Drude model.
    
    Parameters
    ----------
    energy_eV : float or np.ndarray
        Energy values in eV
        
    Returns
    -------
    Material
        Silver material object
    """
    energy_array = np.asarray(energy_eV)
    
    # Simplified Drude model for Ag
    omega_p = 9.0  # Plasma frequency in eV
    gamma = 0.02   # Damping parameter in eV
    
    eps_inf = 1.0  # High-frequency dielectric constant
    eps_ag = eps_inf - (omega_p**2) / (energy_array**2 + 1j * gamma * energy_array)
    
    return Material(
        name="Ag",
        energy_eV=energy_array,
        eps_xx=eps_ag,
        description="Silver with simplified Drude model",
        references="Simplified Drude model - use experimental data for accuracy"
    )

def get_glass(energy_eV: Union[float, np.ndarray], 
              n: float = 1.5,
              material_name: str = "Glass") -> Material:
    """
    Create generic glass material.
    
    Parameters
    ----------
    energy_eV : float or np.ndarray
        Energy values in eV
    n : float, default=1.5
        Refractive index
    material_name : str, default="Glass"
        Name of the glass material
        
    Returns
    -------
    Material
        Glass material object
    """
    energy_array = np.asarray(energy_eV)
    eps_glass = np.full_like(energy_array, n**2, dtype=complex)
    
    return Material(
        name=material_name,
        energy_eV=energy_array,
        eps_xx=eps_glass,
        description=f"Generic glass with n = {n}",
        references="Typical value for common glass"
    )

def get_builtin_material(name: str, 
                        energy_eV: Union[float, np.ndarray],
                        **kwargs) -> Material:
    """
    Get a built-in material by name.
    
    Parameters
    ----------
    name : str
        Name of the material (case-insensitive)
    energy_eV : float or np.ndarray
        Energy values in eV
    **kwargs
        Additional parameters for specific materials
        
    Returns
    -------
    Material
        Requested material object
        
    Raises
    ------
    ValueError
        If material name is not recognized
    """
    name_lower = name.lower()
    
    material_creators = {
        'air': get_air,
        'vacuum': get_vacuum,
        'sio2': get_sio2,
        'silica': get_sio2,
        'si': get_si,
        'silicon': get_si,
        'au': get_au,
        'gold': get_au,
        'ag': get_ag,
        'silver': get_ag,
        'glass': get_glass,
    }
    
    if name_lower in material_creators:
        return material_creators[name_lower](energy_eV, **kwargs)
    else:
        available = ', '.join(material_creators.keys())
        raise ValueError(f"Unknown material '{name}'. Available materials: {available}")

def list_builtin_materials() -> list:
    """
    List all available built-in materials.
    
    Returns
    -------
    list
        List of available material names
    """
    return ['Air', 'Vacuum', 'SiO2', 'Si', 'Au', 'Ag', 'Glass']

def get_material_info(name: str) -> str:
    """
    Get information about a built-in material.
    
    Parameters
    ----------
    name : str
        Name of the material
        
    Returns
    -------
    str
        Information about the material
    """
    info = {
        'air': "Air with refractive index ≈ 1.0",
        'vacuum': "Perfect vacuum with refractive index = 1.0",
        'sio2': "Silicon dioxide (SiO₂) with typical n = 1.46",
        'si': "Silicon with optional dispersion model",
        'au': "Gold with simplified Drude model",
        'ag': "Silver with simplified Drude model", 
        'glass': "Generic glass with configurable refractive index"
    }
    
    name_lower = name.lower()
    if name_lower in info:
        return info[name_lower]
    else:
        return f"No information available for material '{name}'" 