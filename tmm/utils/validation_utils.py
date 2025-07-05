"""
Validation utilities for TMM calculations.

This module provides input validation and error checking functions
to ensure robust operation of TMM calculations.
"""

import numpy as np
from typing import Union, List, Tuple, Any
import warnings

class ValidationError(Exception):
    """Custom exception for validation errors."""
    pass

class ValidationUtils:
    """
    Collection of validation utilities for TMM calculations.
    """
    
    @staticmethod
    def validate_energy(energy: Union[float, np.ndarray], min_val: float = 0.01, max_val: float = 100.0) -> None:
        """
        Validate photon energy values.
        
        Parameters
        ----------
        energy : float or np.ndarray
            Energy in eV
        min_val : float, default=0.01
            Minimum allowed energy in eV
        max_val : float, default=100.0
            Maximum allowed energy in eV
            
        Raises
        ------
        ValidationError
            If energy values are outside allowed range
        """
        energy_arr = np.asarray(energy)
        
        if np.any(energy_arr <= 0):
            raise ValidationError("Energy values must be positive")
        
        if np.any(energy_arr < min_val):
            raise ValidationError(f"Energy values must be >= {min_val} eV")
        
        if np.any(energy_arr > max_val):
            warnings.warn(f"Energy values > {max_val} eV may be outside typical optical range")
    
    @staticmethod
    def validate_thickness(thickness: Union[float, np.ndarray], min_val: float = 0.0, max_val: float = 1e-3) -> None:
        """
        Validate layer thickness values.
        
        Parameters
        ----------
        thickness : float or np.ndarray
            Thickness in meters
        min_val : float, default=0.0
            Minimum allowed thickness in meters
        max_val : float, default=1e-3
            Maximum allowed thickness in meters (1 mm)
            
        Raises
        ------
        ValidationError
            If thickness values are outside allowed range
        """
        thickness_arr = np.asarray(thickness)
        
        if np.any(thickness_arr < min_val):
            raise ValidationError(f"Thickness values must be >= {min_val} m")
        
        if np.any(thickness_arr > max_val):
            warnings.warn(f"Thickness values > {max_val} m may be outside typical film range")
    
    @staticmethod
    def validate_dielectric_function(eps: Union[complex, np.ndarray]) -> None:
        """
        Validate dielectric function values.
        
        Parameters
        ----------
        eps : complex or np.ndarray
            Dielectric function values
            
        Raises
        ------
        ValidationError
            If dielectric function values are unphysical
        """
        eps_arr = np.asarray(eps)
        
        # Check for NaN or infinite values
        if np.any(~np.isfinite(eps_arr)):
            raise ValidationError("Dielectric function contains NaN or infinite values")
        
        # Real part should be positive for most materials
        if np.any(np.real(eps_arr) <= 0):
            warnings.warn("Dielectric function has non-positive real part - may indicate metamaterial")
        
        # Imaginary part should be non-negative for passive materials
        if np.any(np.imag(eps_arr) < 0):
            warnings.warn("Dielectric function has negative imaginary part - may indicate gain medium")
    
    @staticmethod
    def validate_k_vector(kx: Union[float, np.ndarray], ky: Union[float, np.ndarray], 
                         k_max: float = 1e8) -> None:
        """
        Validate k-vector components.
        
        Parameters
        ----------
        kx : float or np.ndarray
            x-component of k-vector in m⁻¹
        ky : float or np.ndarray
            y-component of k-vector in m⁻¹
        k_max : float, default=1e8
            Maximum allowed k-vector magnitude in m⁻¹
            
        Raises
        ------
        ValidationError
            If k-vector components are outside allowed range
        """
        kx_arr = np.asarray(kx)
        ky_arr = np.asarray(ky)
        
        if np.any(~np.isfinite(kx_arr)) or np.any(~np.isfinite(ky_arr)):
            raise ValidationError("k-vector components contain NaN or infinite values")
        
        k_mag = np.sqrt(kx_arr**2 + ky_arr**2)
        if np.any(k_mag > k_max):
            raise ValidationError(f"k-vector magnitude exceeds maximum value {k_max} m⁻¹")
    
    @staticmethod
    def validate_layer_structure(layers: List[Tuple[Any, float]]) -> None:
        """
        Validate layer structure definition.
        
        Parameters
        ----------
        layers : List[Tuple[Any, float]]
            List of (material, thickness) tuples
            
        Raises
        ------
        ValidationError
            If layer structure is invalid
        """
        if len(layers) < 2:
            raise ValidationError("Layer structure must contain at least 2 layers")
        
        # Check first and last layers (typically substrate and superstrate)
        first_thickness = layers[0][1]
        last_thickness = layers[-1][1]
        
        if first_thickness < 1e-3:  # Less than 1 mm
            warnings.warn("First layer thickness < 1 mm - may not represent semi-infinite substrate")
        
        if last_thickness < 1e-3:  # Less than 1 mm
            warnings.warn("Last layer thickness < 1 mm - may not represent semi-infinite superstrate")
        
        # Check intermediate layers
        for i, (material, thickness) in enumerate(layers[1:-1], 1):
            if thickness <= 0:
                raise ValidationError(f"Layer {i} has non-positive thickness")
            
            if thickness > 1e-3:  # Greater than 1 mm
                warnings.warn(f"Layer {i} thickness > 1 mm - unusually thick for thin film")
    
    @staticmethod
    def validate_polarization_parameters(a: float, b: float, delta_phi: float) -> None:
        """
        Validate polarization parameters.
        
        Parameters
        ----------
        a : float
            s-polarization amplitude
        b : float
            p-polarization amplitude
        delta_phi : float
            Phase difference in degrees
            
        Raises
        ------
        ValidationError
            If polarization parameters are invalid
        """
        if a < 0 or b < 0:
            raise ValidationError("Polarization amplitudes must be non-negative")
        
        if a == 0 and b == 0:
            raise ValidationError("At least one polarization amplitude must be non-zero")
        
        if abs(delta_phi) > 360:
            warnings.warn("Phase difference > 360° - will be reduced modulo 360°")
    
    @staticmethod
    def validate_grid_parameters(start: float, stop: float, step: float) -> None:
        """
        Validate grid parameters for calculations.
        
        Parameters
        ----------
        start : float
            Start value
        stop : float
            Stop value
        step : float
            Step size
            
        Raises
        ------
        ValidationError
            If grid parameters are invalid
        """
        if start >= stop:
            raise ValidationError("Start value must be less than stop value")
        
        if step <= 0:
            raise ValidationError("Step size must be positive")
        
        if step > (stop - start):
            raise ValidationError("Step size must be smaller than range")
        
        num_points = int((stop - start) / step) + 1
        if num_points > 10000:
            warnings.warn(f"Grid has {num_points} points - may be computationally expensive")
    
    @staticmethod
    def check_memory_usage(array_shape: Tuple[int, ...], dtype: type = complex) -> None:
        """
        Estimate memory usage for array calculations.
        
        Parameters
        ----------
        array_shape : tuple
            Shape of array to be created
        dtype : type, default=complex
            Data type of array elements
            
        Raises
        ------
        ValidationError
            If estimated memory usage is too high
        """
        num_elements = np.prod(array_shape)
        
        # Estimate bytes per element
        if dtype == complex:
            bytes_per_element = 16  # 8 bytes each for real and imaginary parts
        elif dtype == float:
            bytes_per_element = 8
        else:
            bytes_per_element = 8  # Conservative estimate
        
        estimated_bytes = num_elements * bytes_per_element
        estimated_mb = estimated_bytes / (1024 * 1024)
        
        if estimated_mb > 1000:  # More than 1 GB
            raise ValidationError(f"Estimated memory usage {estimated_mb:.1f} MB exceeds limit")
        elif estimated_mb > 100:  # More than 100 MB
            warnings.warn(f"Estimated memory usage {estimated_mb:.1f} MB is quite large") 