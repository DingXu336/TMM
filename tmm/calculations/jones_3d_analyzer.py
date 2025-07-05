"""
3D Jones Matrix Analyzer - E-kx-ky Analysis

This module implements the full 3D Jones matrix analysis combining energy 
and k-space dependencies. It reproduces and extends the functionality from
the original archive codes with enhanced computational efficiency.

Author: Ding Xu (Physical Chemistry Researcher)
Based on archive/Jones_disperion.py and archive/Jones_Xect.py
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Any
import time
import warnings

from ..core.layer import Layer
from ..materials.material import Material
from ..utils.constants import Constants
from ..utils.physics_utils import PhysicsUtils
from .jones_unified_calculator import JonesUnifiedCalculator, CalculationParams, JonesResults


class Jones3DAnalyzer:
    """
    Specialized 3D Jones Matrix Analyzer for E-kx-ky analysis.
    
    This class implements the complete Jones matrix formalism from the archive
    codes with enhanced computational efficiency and memory management.
    """
    
    def __init__(self, layers: List[Layer]):
        """
        Initialize the 3D analyzer.
        
        Parameters
        ----------
        layers : List[Layer]
            List of layers in the structure (from top to bottom)
        """
        self.layers = layers
        self.validate_layer_structure()
        
        # Physical constants
        self.c = Constants.C
        self.mu0 = Constants.MU0
        self.eps0 = Constants.EPS0
    
    def validate_layer_structure(self) -> None:
        """Validate the layer structure for TMM calculations"""
        if len(self.layers) < 2:
            raise ValueError("Need at least 2 layers (superstrate and substrate)")
        
        # Check for semi-infinite boundaries
        if not self.layers[0].is_semi_infinite:
            warnings.warn("First layer should be semi-infinite (superstrate)")
        if not self.layers[-1].is_semi_infinite:
            warnings.warn("Last layer should be semi-infinite (substrate)")
    
    def get_material_eps(self, layer_index: int, energy_eV: float) -> np.ndarray:
        """
        Get dielectric tensor for a layer at given energy.
        
        Parameters
        ----------
        layer_index : int
            Index of the layer
        energy_eV : float
            Energy in eV
            
        Returns
        -------
        np.ndarray
            3-element dielectric tensor [eps_x, eps_y, eps_z]
        """
        layer = self.layers[layer_index]
        eps_tensor = layer.get_dielectric_tensor(energy_eV)
        
        # Convert to the format expected by the Jones matrix code
        if eps_tensor.ndim == 1:
            return eps_tensor
        else:
            # Extract diagonal elements
            return np.array([eps_tensor[0], eps_tensor[1], eps_tensor[2]])
    
    def local_basis_rotation(self, kx: float, ky: float) -> np.ndarray:
        """
        Calculate local coordinate system rotation matrix.
        
        Parameters
        ----------
        kx, ky : float
            k-vector components in m⁻¹
            
        Returns
        -------
        np.ndarray
            3×3 rotation matrix [s_hat, p_hat, z_hat]
        """
        k_par = np.hypot(kx, ky)
        if k_par == 0:
            return np.eye(3)
        
        # s-polarization unit vector (perpendicular to k_par)
        s_hat = np.array([-ky, kx, 0]) / k_par
        z_hat = np.array([0, 0, 1])
        # p-polarization unit vector  
        p_hat = np.cross(z_hat, s_hat)
        
        return np.column_stack([s_hat, p_hat, z_hat])
    
    def rotate_eps_tensor(self, eps_vec: np.ndarray, kx: float, ky: float) -> np.ndarray:
        """
        Rotate dielectric tensor to local coordinate system.
        
        Parameters
        ----------
        eps_vec : np.ndarray
            Diagonal dielectric tensor [eps_x, eps_y, eps_z]
        kx, ky : float
            k-vector components in m⁻¹
            
        Returns
        -------
        np.ndarray
            3×3 rotated dielectric tensor
        """
        R = self.local_basis_rotation(kx, ky)
        eps_mat = np.diag(eps_vec)
        return R.T @ eps_mat @ R
    
    def interface_jones_matrix(self, eps1: np.ndarray, eps2: np.ndarray, 
                              k_par: float, omega: float) -> np.ndarray:
        """
        Calculate Jones matrix for interface between two media.
        
        Parameters
        ----------
        eps1, eps2 : np.ndarray
            3×3 dielectric tensors of media 1 and 2
        k_par : float
            Parallel component of k-vector in m⁻¹
        omega : float
            Angular frequency in rad/s
            
        Returns
        -------
        np.ndarray
            2×2 Jones reflection matrix
        """
        # Calculate z-components of k-vectors for s and p polarizations
        kz1s = np.sqrt(eps1[1,1] * (omega/self.c)**2 - k_par**2 + 0j)
        kz2s = np.sqrt(eps2[1,1] * (omega/self.c)**2 - k_par**2 + 0j)
        kz1p = np.sqrt(eps1[0,0] * (omega/self.c)**2 - (eps1[0,0]/eps1[2,2]) * k_par**2 + 0j)
        kz2p = np.sqrt(eps2[0,0] * (omega/self.c)**2 - (eps2[0,0]/eps2[2,2]) * k_par**2 + 0j)
        
        # Calculate optical admittances
        Ys1 = kz1s / (self.mu0 * omega)
        Yp1 = self.eps0 * eps1[2,2] * omega / kz1p
        Ys2 = kz2s / (self.mu0 * omega)  
        Yp2 = self.eps0 * eps2[2,2] * omega / kz2p
        
        Y1 = np.diag([Ys1, Yp1])
        Y2 = np.diag([Ys2, Yp2])
        
        # Jones reflection matrix: R = (Y2 - Y1) / (Y2 + Y1)
        return np.linalg.solve(Y2 + Y1, Y2 - Y1)
    
    def gamma_jones_propagation(self, rj: np.ndarray, Gp: np.ndarray, 
                               phi: complex) -> np.ndarray:
        """
        Propagate Jones matrix through layer with phase.
        
        Parameters
        ----------
        rj : np.ndarray
            Interface reflection matrix
        Gp : np.ndarray
            Previous Jones matrix
        phi : complex
            Propagation phase
            
        Returns
        -------
        np.ndarray
            Updated Jones matrix after propagation
        """
        exp2 = np.exp(2j * phi)
        num = rj + Gp * exp2
        den = np.eye(2) + rj @ (Gp * exp2)
        return np.linalg.solve(den, num)
    
    def build_jones_matrix(self, energy_eV: float, kx: float, ky: float) -> Optional[np.ndarray]:
        """
        Build complete Jones matrix for multilayer structure.
        
        Parameters
        ----------
        energy_eV : float
            Energy in eV
        kx, ky : float
            k-vector components in m⁻¹
            
        Returns
        -------
        np.ndarray or None
            2×2 Jones matrix for the complete structure
        """
        # Convert energy to frequency
        omega = energy_eV * Constants.EV_TO_J / Constants.HBAR
        k_par = np.hypot(kx, ky)
        
        # Start from the bottom interface and work upward
        G = None
        
        for iface in reversed(range(len(self.layers) - 1)):
            # Get dielectric tensors for this interface
            eps1 = self.get_material_eps(iface, energy_eV)
            eps2 = self.get_material_eps(iface + 1, energy_eV)
            
            # Rotate to local coordinate system
            eps1_rot = self.rotate_eps_tensor(eps1, kx, ky)
            eps2_rot = self.rotate_eps_tensor(eps2, kx, ky)
            
            # Calculate interface Jones matrix
            rj = self.interface_jones_matrix(eps1_rot, eps2_rot, k_par, omega)
            
            if iface == len(self.layers) - 2:  # Last (bottom) interface
                G = rj
            else:
                # Calculate propagation phase through the layer
                eps_layer = eps2_rot
                kz_p = np.sqrt(eps_layer[0,0] * (omega/self.c)**2 - 
                              (eps_layer[0,0]/eps_layer[2,2]) * k_par**2 + 0j)
                
                # Layer thickness
                layer_thickness = self.layers[iface + 1].thickness
                phi_prop = kz_p * layer_thickness
                
                # Propagate through layer
                G = self.gamma_jones_propagation(rj, G, phi_prop)
        
        return G
    
    def calculate_stokes_parameters(self, Ex: complex, Ey: complex) -> Tuple[float, float, float, float]:
        """
        Calculate Stokes parameters from field components.
        
        Parameters
        ----------
        Ex, Ey : complex
            Electric field components
            
        Returns
        -------
        Tuple[float, float, float, float]
            Stokes parameters (S0, S1, S2, S3)
        """
        S0 = abs(Ex)**2 + abs(Ey)**2
        if S0 == 0:
            return 0.0, 0.0, 0.0, 0.0
        
        S1 = (abs(Ex)**2 - abs(Ey)**2) / S0
        S2 = 2 * np.real(Ex * np.conj(Ey)) / S0
        S3 = 2 * np.imag(Ex * np.conj(Ey)) / S0
        
        return S0, S1, S2, S3
    
    def calculate_polarization_angles(self, S0: float, S1: float, S2: float, S3: float) -> Tuple[float, float]:
        """
        Calculate polarization ellipse angles from Stokes parameters.
        
        Parameters
        ----------
        S0, S1, S2, S3 : float
            Stokes parameters
            
        Returns
        -------
        Tuple[float, float]
            Polarization angles (phi, chi)
        """
        if S0 == 0:
            return 0.0, 0.0
        
        phi = 0.5 * np.arctan2(S2, S1)
        chi = 0.5 * np.arcsin(np.clip(S3 / S0, -1, 1))
        
        return phi, chi
    
    def calculate_single_point(self, energy_eV: float, kx_um: float, ky_um: float,
                             es_amplitude: float, ep_amplitude: float, 
                             phase_diff_deg: float) -> Dict[str, float]:
        """
        Calculate all optical quantities for a single (E, kx, ky) point.
        
        Parameters
        ----------
        energy_eV : float
            Energy in eV
        kx_um, ky_um : float
            k-vector components in μm⁻¹
        es_amplitude, ep_amplitude : float
            Polarization amplitudes
        phase_diff_deg : float
            Phase difference in degrees
            
        Returns
        -------
        Dict[str, float]
            Dictionary with all optical quantities
        """
        # Convert units
        kx = kx_um * 1e6  # μm⁻¹ to m⁻¹
        ky = ky_um * 1e6
        
        # Build Jones matrix
        G = self.build_jones_matrix(energy_eV, kx, ky)
        
        if G is None:
            # Return zeros if calculation failed
            return {
                'Rs': 0.0, 'Rp': 0.0,
                'S0': 0.0, 'S1': 0.0, 'S2': 0.0, 'S3': 0.0,
                'phi': 0.0, 'chi': 0.0
            }
        
        # Input polarization in lab frame
        delta = np.deg2rad(phase_diff_deg)
        Ein_lab = np.array([-ep_amplitude * np.exp(1j * delta), es_amplitude, 0])
        
        # Transform to local coordinates
        R = self.local_basis_rotation(kx, ky)
        Ein_loc = R.T @ Ein_lab
        
        # Apply Jones matrix
        Eout_loc = G @ Ein_loc[:2]
        
        # Transform back to lab frame
        Eref_lab = R @ np.array([Eout_loc[0], Eout_loc[1], 0])
        Ex, Ey = Eref_lab[0], Eref_lab[1]
        
        # Calculate optical quantities
        S0, S1, S2, S3 = self.calculate_stokes_parameters(Ex, Ey)
        phi, chi = self.calculate_polarization_angles(S0, S1, S2, S3)
        
        Rs = abs(G[0, 0])**2  # s-polarized reflectivity
        Rp = abs(G[1, 1])**2  # p-polarized reflectivity
        
        return {
            'Rs': Rs, 'Rp': Rp,
            'S0': S0, 'S1': S1, 'S2': S2, 'S3': S3,
            'phi': phi, 'chi': chi
        }
    
    def calculate_3d_full(self, params: CalculationParams) -> JonesResults:
        """
        Calculate complete 3D E-kx-ky analysis.
        
        Parameters
        ----------
        params : CalculationParams
            Calculation parameters
            
        Returns
        -------
        JonesResults
            Complete 3D analysis results
        """
        start_time = time.time()
        
        # Create coordinate grids
        energy_eV = np.linspace(params.energy_start_eV, params.energy_stop_eV, params.energy_points)
        kx_um = np.linspace(params.kx_start_um, params.kx_stop_um, params.kx_points)
        ky_um = np.linspace(params.ky_start_um, params.ky_stop_um, params.ky_points)
        
        # Initialize result arrays
        shape = (params.energy_points, params.kx_points, params.ky_points)
        Rs = np.zeros(shape, dtype=np.float64)
        Rp = np.zeros(shape, dtype=np.float64)
        S0 = np.zeros(shape, dtype=np.float64)
        S1 = np.zeros(shape, dtype=np.float64)
        S2 = np.zeros(shape, dtype=np.float64)
        S3 = np.zeros(shape, dtype=np.float64)
        phi = np.zeros(shape, dtype=np.float64)
        chi = np.zeros(shape, dtype=np.float64)
        
        total_points = params.energy_points * params.kx_points * params.ky_points
        completed_points = 0
        
        print(f"Starting 3D calculation: {params.energy_points} × {params.kx_points} × {params.ky_points} = {total_points} points")
        
        # Main calculation loop
        for i, E in enumerate(energy_eV):
            print(f"Energy {i+1}/{params.energy_points}: {E:.3f} eV")
            
            for j, kx in enumerate(kx_um):
                for k, ky in enumerate(ky_um):
                    # Calculate for this point
                    result = self.calculate_single_point(
                        E, kx, ky, 
                        params.es_amplitude, params.ep_amplitude, 
                        params.phase_diff_deg
                    )
                    
                    # Store results
                    Rs[i, j, k] = result['Rs']
                    Rp[i, j, k] = result['Rp']
                    S0[i, j, k] = result['S0']
                    S1[i, j, k] = result['S1']
                    S2[i, j, k] = result['S2']
                    S3[i, j, k] = result['S3']
                    phi[i, j, k] = result['phi']
                    chi[i, j, k] = result['chi']
                    
                    completed_points += 1
                    
                    # Progress update
                    if completed_points % 1000 == 0:
                        progress = 100 * completed_points / total_points
                        elapsed = time.time() - start_time
                        eta = elapsed * (total_points - completed_points) / completed_points
                        print(f"  Progress: {progress:.1f}% ({completed_points}/{total_points}) "
                              f"ETA: {eta:.1f}s")
        
        computation_time = time.time() - start_time
        memory_usage = self.estimate_memory_usage(params)
        
        print(f"3D calculation completed in {computation_time:.2f} seconds")
        print(f"Memory usage: {memory_usage:.1f} MB")
        
        return JonesResults(
            energy_eV=energy_eV,
            kx_um=kx_um,
            ky_um=ky_um,
            Rs=Rs, Rp=Rp,
            S0=S0, S1=S1, S2=S2, S3=S3,
            phi=phi, chi=chi,
            calculation_type="3D E-kx-ky full analysis",
            grid_resolution=params.resolution.value,
            computation_time=computation_time,
            memory_usage_mb=memory_usage
        )
    
    def estimate_memory_usage(self, params: CalculationParams) -> float:
        """Estimate memory usage in MB"""
        total_points = params.energy_points * params.kx_points * params.ky_points
        # 8 quantities × 8 bytes per float64 = 64 bytes per point
        memory_mb = total_points * 64 / (1024 * 1024)
        return memory_mb
    
    def get_2d_slice(self, results: JonesResults, energy_eV: float) -> Dict[str, np.ndarray]:
        """
        Extract 2D k-space slice from 3D results at specified energy.
        
        Parameters
        ----------
        results : JonesResults
            3D calculation results
        energy_eV : float
            Energy for the slice
            
        Returns
        -------
        Dict[str, np.ndarray]
            2D slices of all quantities
        """
        # Find closest energy index
        energy_idx = np.argmin(np.abs(results.energy_eV - energy_eV))
        
        return {
            'energy_eV': results.energy_eV[energy_idx],
            'kx_um': results.kx_um,
            'ky_um': results.ky_um,
            'Rs': results.Rs[energy_idx, :, :],
            'Rp': results.Rp[energy_idx, :, :],
            'S0': results.S0[energy_idx, :, :],
            'S1': results.S1[energy_idx, :, :],
            'S2': results.S2[energy_idx, :, :],
            'S3': results.S3[energy_idx, :, :],
            'phi': results.phi[energy_idx, :, :],
            'chi': results.chi[energy_idx, :, :]
        }
    
    def get_1d_dispersion(self, results: JonesResults, k_direction: str, 
                         k_value: float) -> Dict[str, np.ndarray]:
        """
        Extract 1D dispersion from 3D results along specified k-direction.
        
        Parameters
        ----------
        results : JonesResults
            3D calculation results
        k_direction : str
            Direction to extract ('kx' or 'ky')
        k_value : float
            k-value for the dispersion line
            
        Returns
        -------
        Dict[str, np.ndarray]
            1D dispersion curves
        """
        if k_direction == 'kx':
            k_idx = np.argmin(np.abs(results.kx_um - k_value))
            return {
                'energy_eV': results.energy_eV,
                'k_um': results.ky_um,
                'Rs': results.Rs[:, k_idx, :],
                'Rp': results.Rp[:, k_idx, :],
                'S0': results.S0[:, k_idx, :],
                'S1': results.S1[:, k_idx, :],
                'S2': results.S2[:, k_idx, :],
                'S3': results.S3[:, k_idx, :],
                'phi': results.phi[:, k_idx, :],
                'chi': results.chi[:, k_idx, :]
            }
        else:  # k_direction == 'ky'
            k_idx = np.argmin(np.abs(results.ky_um - k_value))
            return {
                'energy_eV': results.energy_eV,
                'k_um': results.kx_um,
                'Rs': results.Rs[:, :, k_idx],
                'Rp': results.Rp[:, :, k_idx],
                'S0': results.S0[:, :, k_idx],
                'S1': results.S1[:, :, k_idx],
                'S2': results.S2[:, :, k_idx],
                'S3': results.S3[:, :, k_idx],
                'phi': results.phi[:, :, k_idx],
                'chi': results.chi[:, :, k_idx]
            } 