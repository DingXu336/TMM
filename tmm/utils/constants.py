"""
Physical constants used in TMM calculations.

This module provides commonly used physical constants in SI units,
following CODATA 2018 recommendations where applicable.
"""

import math

class Constants:
    """
    Collection of physical constants used in TMM calculations.
    
    All constants are in SI units unless otherwise specified.
    """
    
    # Speed of light in vacuum (m/s)
    C = 2.99792458e8
    
    # Vacuum permeability (H/m)
    MU0 = 4e-7 * math.pi
    
    # Vacuum permittivity (F/m)
    EPS0 = 1.0 / (MU0 * C**2)
    
    # Planck constant (J·s)
    H = 6.62607015e-34
    
    # Reduced Planck constant (J·s)
    HBAR = H / (2 * math.pi)
    
    # Elementary charge (C)
    E = 1.602176634e-19
    
    # Electron mass (kg)
    ME = 9.1093837015e-31
    
    # Conversion factors
    EV_TO_J = E  # eV to Joules
    CM_TO_M = 1e-2  # cm to meters
    NM_TO_M = 1e-9  # nm to meters
    UM_TO_M = 1e-6  # μm to meters
    
    # Spectroscopic conversions
    EV_TO_WAVENUMBER = 8065.54429  # eV to cm⁻¹
    WAVENUMBER_TO_EV = 1.0 / EV_TO_WAVENUMBER  # cm⁻¹ to eV
    
    @classmethod
    def energy_to_frequency(cls, energy_eV: float) -> float:
        """Convert energy in eV to frequency in Hz."""
        return energy_eV * cls.EV_TO_J / cls.H
    
    @classmethod
    def energy_to_wavelength(cls, energy_eV: float) -> float:
        """Convert energy in eV to wavelength in meters."""
        frequency = cls.energy_to_frequency(energy_eV)
        return cls.C / frequency
    
    @classmethod
    def energy_to_wavenumber(cls, energy_eV: float) -> float:
        """Convert energy in eV to wavenumber in cm⁻¹."""
        return energy_eV * cls.EV_TO_WAVENUMBER
    
    @classmethod
    def wavenumber_to_energy(cls, wavenumber_cm: float) -> float:
        """Convert wavenumber in cm⁻¹ to energy in eV."""
        return wavenumber_cm * cls.WAVENUMBER_TO_EV
    
    @classmethod
    def wavelength_to_wavenumber(cls, wavelength_um: float) -> float:
        """Convert wavelength in μm to wavenumber in cm⁻¹."""
        return 1e4 / wavelength_um  # 1 μm = 1e-4 cm 