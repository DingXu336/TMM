"""
Layer class for TMM calculations.

This module defines the Layer class that represents individual layers
in a multilayer optical structure.
"""

import numpy as np
from typing import Union, Optional, Dict, Any
from ..materials.material import Material
from ..utils.validation_utils import ValidationUtils

class Layer:
    """
    Represents a single layer in a multilayer optical structure.
    
    A layer consists of a material and a thickness. The layer can be
    finite (with specified thickness) or semi-infinite (substrate/superstrate).
    
    Parameters
    ----------
    material : Material
        The material of the layer
    thickness : float
        Layer thickness in meters (use np.inf for semi-infinite layers)
    name : str, optional
        Name of the layer (defaults to material name)
    """
    
    def __init__(self, 
                 material: Material, 
                 thickness: float,
                 name: Optional[str] = None):
        
        self.material = material
        self.thickness = thickness
        self.name = name or material.name
        
        # Validate inputs
        self._validate_layer()
        
        # Determine if layer is semi-infinite
        self.is_semi_infinite = np.isinf(thickness)
        
        # Store thickness in nanometers for convenience
        self.thickness_nm = thickness * 1e9 if not self.is_semi_infinite else np.inf
    
    def _validate_layer(self) -> None:
        """Validate layer parameters."""
        if not isinstance(self.material, Material):
            raise TypeError("material must be a Material object")
        
        if self.thickness < 0:
            raise ValueError("thickness must be non-negative")
        
        if not np.isinf(self.thickness):
            ValidationUtils.validate_thickness(self.thickness)
    
    def get_dielectric_tensor(self, energy_eV: Union[float, np.ndarray]) -> np.ndarray:
        """
        Get dielectric tensor at specified energies.
        
        Parameters
        ----------
        energy_eV : float or np.ndarray
            Energy values in eV
            
        Returns
        -------
        np.ndarray
            Dielectric tensor with shape (..., 3)
        """
        return self.material.get_dielectric_tensor(energy_eV)
    
    def get_optical_constants(self, energy_eV: Union[float, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Get optical constants at specified energies.
        
        Parameters
        ----------
        energy_eV : float or np.ndarray
            Energy values in eV
            
        Returns
        -------
        Dict[str, np.ndarray]
            Dictionary with optical constants
        """
        return self.material.get_optical_constants(energy_eV)
    
    def get_propagation_phase(self, energy_eV: Union[float, np.ndarray], 
                             kz: Union[complex, np.ndarray]) -> Union[complex, np.ndarray]:
        """
        Calculate propagation phase through the layer.
        
        Parameters
        ----------
        energy_eV : float or np.ndarray
            Energy values in eV
        kz : complex or np.ndarray
            z-component of k-vector
            
        Returns
        -------
        complex or np.ndarray
            Propagation phase factor exp(i*kz*d)
        """
        if self.is_semi_infinite:
            # Semi-infinite layers don't accumulate phase
            return np.ones_like(kz, dtype=complex)
        
        return np.exp(1j * kz * self.thickness)
    
    def get_optical_thickness(self, energy_eV: Union[float, np.ndarray],
                             angle_deg: float = 0.0) -> Union[float, np.ndarray]:
        """
        Calculate optical thickness of the layer.
        
        Parameters
        ----------
        energy_eV : float or np.ndarray
            Energy values in eV
        angle_deg : float, default=0.0
            Angle of incidence in degrees
            
        Returns
        -------
        float or np.ndarray
            Optical thickness in units of wavelength
        """
        if self.is_semi_infinite:
            return np.inf
        
        # Get refractive index
        n_k = self.material.get_refractive_index(energy_eV)
        n_eff = n_k['n_x']  # Use x-component for simplicity
        
        # Calculate wavelength in material
        from ..utils.constants import Constants
        wavelength_vacuum = Constants.energy_to_wavelength(energy_eV)
        
        # Account for angle of incidence
        cos_theta = np.cos(np.deg2rad(angle_deg))
        
        optical_thickness = self.thickness * n_eff / (wavelength_vacuum * cos_theta)
        
        return optical_thickness
    
    def is_absorbing(self, energy_eV: Union[float, np.ndarray], 
                    threshold: float = 1e-6) -> Union[bool, np.ndarray]:
        """
        Check if layer is absorbing at given energies.
        
        Parameters
        ----------
        energy_eV : float or np.ndarray
            Energy values in eV
        threshold : float, default=1e-6
            Threshold for considering material as absorbing
            
        Returns
        -------
        bool or np.ndarray
            True if material is absorbing (k > threshold)
        """
        n_k = self.material.get_refractive_index(energy_eV)
        k = n_k['k_x']  # Use x-component for simplicity
        
        return k > threshold
    
    def copy(self) -> 'Layer':
        """
        Create a copy of the layer.
        
        Returns
        -------
        Layer
            Copy of the layer
        """
        return Layer(
            material=self.material,  # Materials are immutable, so we can share
            thickness=self.thickness,
            name=self.name
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert layer to dictionary for serialization.
        
        Returns
        -------
        Dict[str, Any]
            Layer data as dictionary
        """
        return {
            'material': self.material.to_dict(),
            'thickness': self.thickness,
            'name': self.name,
            'is_semi_infinite': self.is_semi_infinite
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Layer':
        """
        Create layer from dictionary.
        
        Parameters
        ----------
        data : Dict[str, Any]
            Layer data dictionary
            
        Returns
        -------
        Layer
            Layer instance
        """
        material = Material.from_dict(data['material'])
        return cls(
            material=material,
            thickness=data['thickness'],
            name=data.get('name')
        )
    
    def __str__(self) -> str:
        """String representation of the layer."""
        if self.is_semi_infinite:
            thickness_str = "semi-infinite"
        else:
            thickness_str = f"{self.thickness_nm:.1f} nm"
        
        return f"Layer(name='{self.name}', material='{self.material.name}', thickness={thickness_str})"
    
    def __repr__(self) -> str:
        """Detailed string representation."""
        return self.__str__()
    
    def __eq__(self, other) -> bool:
        """Check equality with another layer."""
        if not isinstance(other, Layer):
            return False
        
        return (self.material == other.material and
                self.thickness == other.thickness and
                self.name == other.name)

def create_substrate(material: Material, name: Optional[str] = None) -> Layer:
    """
    Create a semi-infinite substrate layer.
    
    Parameters
    ----------
    material : Material
        The substrate material
    name : str, optional
        Name of the substrate layer
        
    Returns
    -------
    Layer
        Semi-infinite substrate layer
    """
    return Layer(
        material=material,
        thickness=np.inf,
        name=name or f"{material.name}_substrate"
    )

def create_superstrate(material: Material, name: Optional[str] = None) -> Layer:
    """
    Create a semi-infinite superstrate layer.
    
    Parameters
    ----------
    material : Material
        The superstrate material
    name : str, optional
        Name of the superstrate layer
        
    Returns
    -------
    Layer
        Semi-infinite superstrate layer
    """
    return Layer(
        material=material,
        thickness=np.inf,
        name=name or f"{material.name}_superstrate"
    )

def create_film(material: Material, thickness_nm: float, 
               name: Optional[str] = None) -> Layer:
    """
    Create a finite film layer.
    
    Parameters
    ----------
    material : Material
        The film material
    thickness_nm : float
        Film thickness in nanometers
    name : str, optional
        Name of the film layer
        
    Returns
    -------
    Layer
        Finite film layer
    """
    thickness_m = thickness_nm * 1e-9
    return Layer(
        material=material,
        thickness=thickness_m,
        name=name or f"{material.name}_film"
    ) 