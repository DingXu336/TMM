"""
Material class for optical materials in TMM calculations.

This module defines the Material class that represents optical materials
with their dielectric properties as functions of energy/wavelength.
"""

import numpy as np
from typing import Union, Optional, Dict, Any
from ..utils.constants import Constants
from ..utils.physics_utils import PhysicsUtils
from ..utils.validation_utils import ValidationUtils

class Material:
    """
    Represents an optical material with dielectric properties.
    
    This class stores and manages the dielectric function (permittivity)
    of a material as a function of energy or wavelength. It supports
    both isotropic and anisotropic materials.
    
    Parameters
    ----------
    name : str
        Name of the material
    energy_eV : np.ndarray
        Energy grid in eV
    eps_xx : np.ndarray
        xx-component of dielectric tensor
    eps_yy : np.ndarray, optional
        yy-component of dielectric tensor (defaults to eps_xx)
    eps_zz : np.ndarray, optional
        zz-component of dielectric tensor (defaults to eps_xx)
    description : str, optional
        Description of the material
    references : str, optional
        Literature references for the material data
    """
    
    def __init__(self, 
                 name: str,
                 energy_eV: Union[np.ndarray, float],
                 eps_xx: Union[np.ndarray, complex],
                 eps_yy: Optional[Union[np.ndarray, complex]] = None,
                 eps_zz: Optional[Union[np.ndarray, complex]] = None,
                 description: Optional[str] = None,
                 references: Optional[str] = None):
        
        self.name = name
        self.description = description or f"Material: {name}"
        self.references = references or "No references provided"
        
        # Convert to arrays
        self.energy_eV = np.asarray(energy_eV)
        self.eps_xx = np.asarray(eps_xx)
        
        # Handle isotropic case
        if eps_yy is None:
            self.eps_yy = self.eps_xx.copy()
        else:
            self.eps_yy = np.asarray(eps_yy)
            
        if eps_zz is None:
            self.eps_zz = self.eps_xx.copy()
        else:
            self.eps_zz = np.asarray(eps_zz)
        
        # Validate inputs
        self._validate_data()
        
        # Check if material is anisotropic
        self.is_anisotropic = not (np.allclose(self.eps_xx, self.eps_yy) and 
                                  np.allclose(self.eps_xx, self.eps_zz))
    
    def _validate_data(self) -> None:
        """Validate material data for consistency."""
        # Validate energy grid
        ValidationUtils.validate_energy(self.energy_eV)
        
        # Validate dielectric functions
        ValidationUtils.validate_dielectric_function(self.eps_xx)
        ValidationUtils.validate_dielectric_function(self.eps_yy)
        ValidationUtils.validate_dielectric_function(self.eps_zz)
        
        # Check array shapes
        if self.energy_eV.ndim > 1:
            raise ValueError("Energy array must be 1D")
        
        # For arrays, check shape consistency
        if self.energy_eV.size > 1:
            if self.eps_xx.shape != self.energy_eV.shape:
                raise ValueError("eps_xx shape must match energy array shape")
            if self.eps_yy.shape != self.energy_eV.shape:
                raise ValueError("eps_yy shape must match energy array shape")
            if self.eps_zz.shape != self.energy_eV.shape:
                raise ValueError("eps_zz shape must match energy array shape")
    
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
            Dielectric tensor with shape (..., 3, 3) for anisotropic materials
            or (..., 3) for diagonal tensors
        """
        energy_array = np.asarray(energy_eV)
        
        # Interpolate dielectric function values
        eps_xx_interp = self._interpolate_dielectric(energy_array, self.eps_xx)
        eps_yy_interp = self._interpolate_dielectric(energy_array, self.eps_yy)
        eps_zz_interp = self._interpolate_dielectric(energy_array, self.eps_zz)
        
        # Return diagonal components for simplicity
        # Full tensor support can be added later if needed
        return np.stack([eps_xx_interp, eps_yy_interp, eps_zz_interp], axis=-1)
    
    def _interpolate_dielectric(self, energy_target: np.ndarray, eps_source: np.ndarray) -> np.ndarray:
        """
        Interpolate dielectric function to target energy grid.
        
        Parameters
        ----------
        energy_target : np.ndarray
            Target energy values
        eps_source : np.ndarray
            Source dielectric function values
            
        Returns
        -------
        np.ndarray
            Interpolated dielectric function values
        """
        if self.energy_eV.size == 1:
            # Single energy point - return constant values
            return np.full_like(energy_target, eps_source.item(), dtype=complex)
        
        # Interpolate real and imaginary parts separately
        eps_real = np.interp(energy_target, self.energy_eV, eps_source.real)
        eps_imag = np.interp(energy_target, self.energy_eV, eps_source.imag)
        
        return eps_real + 1j * eps_imag
    
    def get_refractive_index(self, energy_eV: Union[float, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Calculate refractive index from dielectric function.
        
        Parameters
        ----------
        energy_eV : float or np.ndarray
            Energy values in eV
            
        Returns
        -------
        Dict[str, np.ndarray]
            Dictionary with 'n' (real part) and 'k' (imaginary part) for each axis
        """
        eps_tensor = self.get_dielectric_tensor(energy_eV)
        
        result = {}
        for i, axis in enumerate(['x', 'y', 'z']):
            eps = eps_tensor[..., i]
            n_complex = np.sqrt(eps)
            result[f'n_{axis}'] = n_complex.real
            result[f'k_{axis}'] = n_complex.imag
        
        return result
    
    def get_optical_constants(self, energy_eV: Union[float, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Get comprehensive optical constants.
        
        Parameters
        ----------
        energy_eV : float or np.ndarray
            Energy values in eV
            
        Returns
        -------
        Dict[str, np.ndarray]
            Dictionary with optical constants: eps, n, k, alpha, etc.
        """
        energy_array = np.asarray(energy_eV)
        eps_tensor = self.get_dielectric_tensor(energy_array)
        n_k = self.get_refractive_index(energy_array)
        
        # Calculate absorption coefficient (m⁻¹)
        wavelength_m = Constants.energy_to_wavelength(energy_array)
        alpha = {}
        for axis in ['x', 'y', 'z']:
            alpha[f'alpha_{axis}'] = 4 * np.pi * n_k[f'k_{axis}'] / wavelength_m
        
        # Combine all results
        result = {
            'energy_eV': energy_array,
            'wavelength_m': wavelength_m,
            'eps_xx': eps_tensor[..., 0],
            'eps_yy': eps_tensor[..., 1],
            'eps_zz': eps_tensor[..., 2],
            **n_k,
            **alpha
        }
        
        return result
    
    def plot_optical_properties(self, energy_range: Optional[tuple] = None, 
                               save_path: Optional[str] = None) -> None:
        """
        Plot optical properties of the material.
        
        Parameters
        ----------
        energy_range : tuple, optional
            (min_energy, max_energy) in eV for plotting
        save_path : str, optional
            Path to save the plot
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib required for plotting")
        
        if energy_range is None:
            energy_plot = self.energy_eV
        else:
            energy_plot = np.linspace(energy_range[0], energy_range[1], 1000)
        
        optical_constants = self.get_optical_constants(energy_plot)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Optical Properties of {self.name}', fontsize=16)
        
        # Dielectric function
        ax = axes[0, 0]
        ax.plot(energy_plot, optical_constants['eps_xx'].real, 'b-', label='Re(εₓₓ)')
        ax.plot(energy_plot, optical_constants['eps_xx'].imag, 'r-', label='Im(εₓₓ)')
        if self.is_anisotropic:
            ax.plot(energy_plot, optical_constants['eps_yy'].real, 'b--', label='Re(εᵧᵧ)')
            ax.plot(energy_plot, optical_constants['eps_zz'].real, 'b:', label='Re(εᵨᵨ)')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Dielectric function')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Refractive index
        ax = axes[0, 1]
        ax.plot(energy_plot, optical_constants['n_x'], 'b-', label='n')
        ax.plot(energy_plot, optical_constants['k_x'], 'r-', label='k')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Refractive index')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Absorption coefficient
        ax = axes[1, 0]
        ax.semilogy(energy_plot, optical_constants['alpha_x'], 'g-', label='α')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Absorption coefficient (m⁻¹)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Wavelength dependence
        ax = axes[1, 1]
        wavelength_nm = optical_constants['wavelength_m'] * 1e9
        ax.plot(wavelength_nm, optical_constants['n_x'], 'b-', label='n')
        ax.plot(wavelength_nm, optical_constants['k_x'], 'r-', label='k')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Refractive index')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        else:
            plt.show()
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert material to dictionary for serialization.
        
        Returns
        -------
        Dict[str, Any]
            Material data as dictionary
        """
        return {
            'name': self.name,
            'description': self.description,
            'references': self.references,
            'energy_eV': self.energy_eV.tolist(),
            'eps_xx': [complex(x) for x in self.eps_xx],
            'eps_yy': [complex(x) for x in self.eps_yy],
            'eps_zz': [complex(x) for x in self.eps_zz],
            'is_anisotropic': self.is_anisotropic
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Material':
        """
        Create material from dictionary.
        
        Parameters
        ----------
        data : Dict[str, Any]
            Material data dictionary
            
        Returns
        -------
        Material
            Material instance
        """
        return cls(
            name=data['name'],
            energy_eV=np.array(data['energy_eV']),
            eps_xx=np.array(data['eps_xx']),
            eps_yy=np.array(data['eps_yy']),
            eps_zz=np.array(data['eps_zz']),
            description=data.get('description'),
            references=data.get('references')
        )
    
    def __str__(self) -> str:
        """String representation of the material."""
        aniso_str = "anisotropic" if self.is_anisotropic else "isotropic"
        return f"Material(name='{self.name}', {aniso_str}, {len(self.energy_eV)} energy points)"
    
    def __repr__(self) -> str:
        """Detailed string representation."""
        return (f"Material(name='{self.name}', "
                f"energy_range=({self.energy_eV.min():.2f}, {self.energy_eV.max():.2f}) eV, "
                f"anisotropic={self.is_anisotropic})") 