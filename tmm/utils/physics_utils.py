"""
Physics utilities for TMM calculations.

This module provides common physics calculations, coordinate transformations,
and mathematical operations used in Transfer Matrix Method calculations.
"""

import numpy as np
from typing import Tuple, Union, Optional
from .constants import Constants

class PhysicsUtils:
    """
    Collection of physics utility functions for TMM calculations.
    """
    
    @staticmethod
    def local_basis_rotation(kx: float, ky: float) -> np.ndarray:
        """
        Compute local basis rotation matrix for given k-vector components.
        
        This function creates a rotation matrix that transforms from the laboratory
        frame to the local coordinate system where the incident k-vector lies in
        the x-z plane.
        
        Parameters
        ----------
        kx : float
            x-component of k-vector in SI units (m⁻¹)
        ky : float
            y-component of k-vector in SI units (m⁻¹)
            
        Returns
        -------
        np.ndarray
            3x3 rotation matrix [s_hat, p_hat, z_hat]
            
        Notes
        -----
        The local basis is defined as:
        - s_hat: perpendicular to incidence plane (TE direction)
        - p_hat: in incidence plane, perpendicular to k (TM direction)  
        - z_hat: surface normal (growth direction)
        """
        k_par = np.sqrt(kx**2 + ky**2)
        
        if k_par == 0:
            # Normal incidence - no rotation needed
            return np.eye(3)
        
        # s-polarization unit vector (TE)
        s_hat = np.array([-ky, kx, 0]) / k_par
        
        # z-direction (surface normal)
        z_hat = np.array([0, 0, 1])
        
        # p-polarization unit vector (TM)
        p_hat = np.cross(z_hat, s_hat)
        
        return np.column_stack([s_hat, p_hat, z_hat])
    
    @staticmethod
    def rotate_dielectric_tensor(eps_diag: np.ndarray, kx: float, ky: float) -> np.ndarray:
        """
        Rotate diagonal dielectric tensor to local coordinate system.
        
        Parameters
        ----------
        eps_diag : np.ndarray
            Diagonal dielectric tensor components [eps_x, eps_y, eps_z]
        kx : float
            x-component of k-vector in SI units (m⁻¹)
        ky : float
            y-component of k-vector in SI units (m⁻¹)
            
        Returns
        -------
        np.ndarray
            3x3 rotated dielectric tensor
        """
        R = PhysicsUtils.local_basis_rotation(kx, ky)
        eps_matrix = np.diag(eps_diag)
        return R.T @ eps_matrix @ R
    
    @staticmethod
    def compute_kz_components(eps_tensor: np.ndarray, k_parallel: float, omega: float) -> Tuple[complex, complex]:
        """
        Compute z-components of k-vector for s and p polarizations.
        
        Parameters
        ----------
        eps_tensor : np.ndarray
            3x3 dielectric tensor in local coordinates
        k_parallel : float
            Parallel component of k-vector (m⁻¹)
        omega : float
            Angular frequency (rad/s)
            
        Returns
        -------
        Tuple[complex, complex]
            (kz_s, kz_p) - z-components for s and p polarizations
        """
        k0 = omega / Constants.C
        
        # s-polarization (TE): kz² = εyy * k0² - kpar²
        kz_s = np.sqrt(eps_tensor[1, 1] * k0**2 - k_parallel**2 + 0j)
        
        # p-polarization (TM): kz² = εxx * k0² - (εxx/εzz) * kpar²
        kz_p = np.sqrt(eps_tensor[0, 0] * k0**2 - (eps_tensor[0, 0]/eps_tensor[2, 2]) * k_parallel**2 + 0j)
        
        return kz_s, kz_p
    
    @staticmethod
    def compute_admittances(eps_tensor: np.ndarray, kz_s: complex, kz_p: complex, omega: float) -> Tuple[complex, complex]:
        """
        Compute optical admittances for s and p polarizations.
        
        Parameters
        ----------
        eps_tensor : np.ndarray
            3x3 dielectric tensor in local coordinates
        kz_s : complex
            z-component of k-vector for s-polarization
        kz_p : complex
            z-component of k-vector for p-polarization
        omega : float
            Angular frequency (rad/s)
            
        Returns
        -------
        Tuple[complex, complex]
            (Y_s, Y_p) - admittances for s and p polarizations
        """
        # s-polarization admittance
        Y_s = kz_s / (Constants.MU0 * omega)
        
        # p-polarization admittance
        Y_p = Constants.EPS0 * eps_tensor[2, 2] * omega / kz_p
        
        return Y_s, Y_p
    
    @staticmethod
    def stokes_parameters(Ex: complex, Ey: complex) -> Tuple[float, float, float, float]:
        """
        Calculate Stokes parameters from electric field components.
        
        Parameters
        ----------
        Ex : complex
            x-component of electric field
        Ey : complex
            y-component of electric field
            
        Returns
        -------
        Tuple[float, float, float, float]
            (S0, S1, S2, S3) - Stokes parameters
        """
        I_x = np.abs(Ex)**2
        I_y = np.abs(Ey)**2
        
        S0 = I_x + I_y  # Total intensity
        S1 = I_x - I_y  # Linear polarization (horizontal - vertical)
        S2 = 2 * np.real(Ex * np.conj(Ey))  # Linear polarization (45° - 135°)
        S3 = 2 * np.imag(Ex * np.conj(Ey))  # Circular polarization
        
        return S0, S1, S2, S3
    
    @staticmethod
    def polarization_ellipse_parameters(S0: float, S1: float, S2: float, S3: float) -> Tuple[float, float]:
        """
        Calculate polarization ellipse parameters from Stokes parameters.
        
        Parameters
        ----------
        S0 : float
            Total intensity
        S1 : float
            Linear polarization parameter
        S2 : float
            Linear polarization parameter
        S3 : float
            Circular polarization parameter
            
        Returns
        -------
        Tuple[float, float]
            (phi, chi) - ellipse orientation angle and ellipticity angle
        """
        # Avoid division by zero
        if S0 == 0:
            return 0.0, 0.0
        
        # Orientation angle of major axis
        phi = 0.5 * np.arctan2(S2, S1)
        
        # Ellipticity angle
        chi = 0.5 * np.arcsin(np.clip(S3 / S0, -1, 1))
        
        return phi, chi
    
    @staticmethod
    def energy_to_angular_frequency(energy_eV: float) -> float:
        """Convert photon energy in eV to angular frequency in rad/s."""
        return energy_eV * Constants.EV_TO_J / Constants.HBAR
    
    @staticmethod
    def refractive_index_to_dielectric(n: Union[float, complex], k: Optional[float] = None) -> complex:
        """
        Convert refractive index to dielectric function.
        
        Parameters
        ----------
        n : float or complex
            Refractive index (real part)
        k : float, optional
            Extinction coefficient (imaginary part)
            
        Returns
        -------
        complex
            Complex dielectric function ε = (n + ik)²
        """
        if k is not None:
            n_complex = n + 1j * k
        else:
            n_complex = n
            
        return n_complex**2 