"""
Unified Jones Matrix Calculator

This module provides a unified interface for Jones matrix calculations combining:
- 3D E-kx-ky analysis (full parameter space)
- 2D kx-ky cross-sections at fixed energy
- 1D dispersion analysis along k-directions
- Adaptive grid resolution for computational efficiency

Author: Ding Xu (Physical Chemistry Researcher)
"""

import numpy as np
from typing import Union, Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from enum import Enum
import warnings

from ..core.layer import Layer
from ..materials.material import Material
from ..utils.constants import Constants
from ..utils.physics_utils import PhysicsUtils
from ..utils.validation_utils import ValidationUtils

class GridResolution(Enum):
    """Grid resolution levels for adaptive calculation"""
    COARSE = "coarse"
    MEDIUM = "medium"
    FINE = "fine"
    CUSTOM = "custom"

@dataclass
class CalculationParams:
    """Parameters for Jones matrix calculations"""
    # Energy parameters
    energy_start_eV: float = 1.0
    energy_stop_eV: float = 3.0
    energy_points: int = 101
    
    # k-space parameters
    kx_start_um: float = 0.0
    kx_stop_um: float = 25.0
    kx_points: int = 101
    ky_start_um: float = 0.0
    ky_stop_um: float = 25.0
    ky_points: int = 101
    
    # Polarization parameters
    es_amplitude: float = 1.0  # s-polarization amplitude
    ep_amplitude: float = 0.0  # p-polarization amplitude
    phase_diff_deg: float = 0.0  # Phase difference between s and p
    
    # Grid resolution
    resolution: GridResolution = GridResolution.MEDIUM
    
    # Memory management
    max_memory_gb: float = 8.0
    use_chunking: bool = True
    chunk_size: int = 1000

@dataclass
class JonesResults:
    """Results from Jones matrix calculations"""
    # Coordinate grids
    energy_eV: np.ndarray
    kx_um: np.ndarray
    ky_um: np.ndarray
    
    # Reflectivity
    Rs: np.ndarray  # s-polarized reflectivity
    Rp: np.ndarray  # p-polarized reflectivity
    
    # Stokes parameters
    S0: np.ndarray  # Total intensity
    S1: np.ndarray  # Linear H-V polarization
    S2: np.ndarray  # Linear ±45° polarization
    S3: np.ndarray  # Circular polarization
    
    # Polarization ellipse parameters
    phi: np.ndarray  # Azimuth angle
    chi: np.ndarray  # Ellipticity angle
    
    # Metadata
    calculation_type: str
    grid_resolution: str
    computation_time: float
    memory_usage_mb: float

class JonesUnifiedCalculator:
    """
    Unified Jones Matrix Calculator for multilayer optical structures.
    
    This class provides comprehensive Jones matrix analysis capabilities:
    - 3D E-kx-ky mapping with adaptive resolution
    - 2D kx-ky cross-sections at fixed energy
    - 1D dispersion analysis along k-directions
    - Smart memory management for large datasets
    """
    
    def __init__(self, layers: List[Layer]):
        """
        Initialize the calculator with layer structure.
        
        Parameters
        ----------
        layers : List[Layer]
            List of layers in the structure (from top to bottom)
        """
        self.layers = layers
        self.validate_layers()
        
        # Default calculation parameters
        self.default_params = CalculationParams()
        
        # Resolution presets
        self.resolution_presets = {
            GridResolution.COARSE: {
                'energy_factor': 0.5,
                'k_factor': 0.5,
                'description': 'Fast preview calculations'
            },
            GridResolution.MEDIUM: {
                'energy_factor': 1.0,
                'k_factor': 1.0,
                'description': 'Balanced quality and speed'
            },
            GridResolution.FINE: {
                'energy_factor': 2.0,
                'k_factor': 2.0,
                'description': 'High-quality final results'
            }
        }
    
    def validate_layers(self) -> None:
        """Validate layer structure"""
        if not self.layers:
            raise ValueError("Layer structure cannot be empty")
        
        if len(self.layers) < 2:
            raise ValueError("Need at least 2 layers (superstrate and substrate)")
        
        # Check for semi-infinite layers at ends
        if not self.layers[0].is_semi_infinite:
            warnings.warn("First layer should be semi-infinite (superstrate)")
        if not self.layers[-1].is_semi_infinite:
            warnings.warn("Last layer should be semi-infinite (substrate)")
    
    def estimate_memory_usage(self, params: CalculationParams) -> float:
        """
        Estimate memory usage for given parameters.
        
        Parameters
        ----------
        params : CalculationParams
            Calculation parameters
            
        Returns
        -------
        float
            Estimated memory usage in MB
        """
        # Each complex number = 16 bytes
        # 8 output quantities × 16 bytes = 128 bytes per k-point per energy
        total_points = params.energy_points * params.kx_points * params.ky_points
        memory_mb = total_points * 128 / (1024 * 1024)
        
        # Add overhead (intermediate calculations, Jones matrices)
        memory_mb *= 3.0
        
        return memory_mb
    
    def adapt_resolution(self, params: CalculationParams) -> CalculationParams:
        """
        Adapt grid resolution based on memory constraints.
        
        Parameters
        ----------
        params : CalculationParams
            Input parameters
            
        Returns
        -------
        CalculationParams
            Adapted parameters
        """
        if params.resolution == GridResolution.CUSTOM:
            return params
        
        # Get resolution factors
        preset = self.resolution_presets[params.resolution]
        
        # Create adapted parameters
        adapted_params = CalculationParams(
            energy_start_eV=params.energy_start_eV,
            energy_stop_eV=params.energy_stop_eV,
            energy_points=max(1, int(params.energy_points * preset['energy_factor'])),
            
            kx_start_um=params.kx_start_um,
            kx_stop_um=params.kx_stop_um,
            kx_points=max(1, int(params.kx_points * preset['k_factor'])),
            
            ky_start_um=params.ky_start_um,
            ky_stop_um=params.ky_stop_um,
            ky_points=max(1, int(params.ky_points * preset['k_factor'])),
            
            es_amplitude=params.es_amplitude,
            ep_amplitude=params.ep_amplitude,
            phase_diff_deg=params.phase_diff_deg,
            
            resolution=params.resolution,
            max_memory_gb=params.max_memory_gb,
            use_chunking=params.use_chunking,
            chunk_size=params.chunk_size
        )
        
        # Check memory constraints
        estimated_memory = self.estimate_memory_usage(adapted_params)
        if estimated_memory > params.max_memory_gb * 1024:
            warnings.warn(f"Estimated memory usage ({estimated_memory:.1f} MB) exceeds limit. "
                         f"Consider using coarser resolution or chunking.")
        
        return adapted_params
    
    def calculate_3d_mapping(self, params: Optional[CalculationParams] = None) -> JonesResults:
        """
        Calculate complete 3D E-kx-ky mapping.
        
        Parameters
        ----------
        params : CalculationParams, optional
            Calculation parameters (uses defaults if None)
            
        Returns
        -------
        JonesResults
            Complete 3D analysis results
        """
        if params is None:
            params = self.default_params
        
        # Adapt resolution
        adapted_params = self.adapt_resolution(params)
        
        print(f"Starting 3D E-kx-ky calculation with {adapted_params.resolution.value} resolution")
        print(f"Grid size: {adapted_params.energy_points} × {adapted_params.kx_points} × {adapted_params.ky_points}")
        
        # Create coordinate grids
        energy_eV = np.linspace(adapted_params.energy_start_eV, adapted_params.energy_stop_eV, adapted_params.energy_points)
        kx_um = np.linspace(adapted_params.kx_start_um, adapted_params.kx_stop_um, adapted_params.kx_points)
        ky_um = np.linspace(adapted_params.ky_start_um, adapted_params.ky_stop_um, adapted_params.ky_points)
        
        # Initialize result arrays
        shape = (adapted_params.energy_points, adapted_params.kx_points, adapted_params.ky_points)
        Rs = np.zeros(shape, dtype=np.float64)
        Rp = np.zeros(shape, dtype=np.float64)
        S0 = np.zeros(shape, dtype=np.float64)
        S1 = np.zeros(shape, dtype=np.float64)
        S2 = np.zeros(shape, dtype=np.float64)
        S3 = np.zeros(shape, dtype=np.float64)
        phi = np.zeros(shape, dtype=np.float64)
        chi = np.zeros(shape, dtype=np.float64)
        
        # Calculate using chunking if enabled
        if adapted_params.use_chunking:
            return self._calculate_3d_chunked(adapted_params, energy_eV, kx_um, ky_um,
                                            Rs, Rp, S0, S1, S2, S3, phi, chi)
        else:
            return self._calculate_3d_direct(adapted_params, energy_eV, kx_um, ky_um,
                                           Rs, Rp, S0, S1, S2, S3, phi, chi)
    
    def _calculate_3d_direct(self, params: CalculationParams, energy_eV: np.ndarray,
                           kx_um: np.ndarray, ky_um: np.ndarray,
                           Rs: np.ndarray, Rp: np.ndarray, S0: np.ndarray, 
                           S1: np.ndarray, S2: np.ndarray, S3: np.ndarray,
                           phi: np.ndarray, chi: np.ndarray) -> JonesResults:
        """Direct 3D calculation without chunking"""
        import time
        start_time = time.time()
        
        total_points = len(energy_eV) * len(kx_um) * len(ky_um)
        completed_points = 0
        
        for i, E in enumerate(energy_eV):
            print(f"Energy progress: {i+1}/{len(energy_eV)} ({E:.2f} eV)")
            
            for j, kx in enumerate(kx_um):
                for k, ky in enumerate(ky_um):
                    # Calculate Jones matrix for this (E, kx, ky) point
                    results = self._calculate_single_point(E, kx, ky, params)
                    
                    # Store results
                    Rs[i, j, k] = results['Rs']
                    Rp[i, j, k] = results['Rp']
                    S0[i, j, k] = results['S0']
                    S1[i, j, k] = results['S1']
                    S2[i, j, k] = results['S2']
                    S3[i, j, k] = results['S3']
                    phi[i, j, k] = results['phi']
                    chi[i, j, k] = results['chi']
                    
                    completed_points += 1
                    
                    # Progress update
                    if completed_points % 1000 == 0:
                        progress = 100 * completed_points / total_points
                        print(f"  Progress: {progress:.1f}% ({completed_points}/{total_points})")
        
        computation_time = time.time() - start_time
        memory_usage = self.estimate_memory_usage(params)
        
        return JonesResults(
            energy_eV=energy_eV,
            kx_um=kx_um,
            ky_um=ky_um,
            Rs=Rs, Rp=Rp,
            S0=S0, S1=S1, S2=S2, S3=S3,
            phi=phi, chi=chi,
            calculation_type="3D E-kx-ky mapping",
            grid_resolution=params.resolution.value,
            computation_time=computation_time,
            memory_usage_mb=memory_usage
        )
    
    def _calculate_3d_chunked(self, params: CalculationParams, energy_eV: np.ndarray,
                            kx_um: np.ndarray, ky_um: np.ndarray,
                            Rs: np.ndarray, Rp: np.ndarray, S0: np.ndarray, 
                            S1: np.ndarray, S2: np.ndarray, S3: np.ndarray,
                            phi: np.ndarray, chi: np.ndarray) -> JonesResults:
        """Chunked 3D calculation for memory management"""
        import time
        start_time = time.time()
        
        print("Using chunked calculation for memory efficiency")
        
        # Process in chunks by energy
        chunk_size = min(params.chunk_size, len(energy_eV))
        
        for chunk_start in range(0, len(energy_eV), chunk_size):
            chunk_end = min(chunk_start + chunk_size, len(energy_eV))
            energy_chunk = energy_eV[chunk_start:chunk_end]
            
            print(f"Processing energy chunk {chunk_start+1}-{chunk_end}/{len(energy_eV)}")
            
            # Calculate chunk
            for i, E in enumerate(energy_chunk):
                global_i = chunk_start + i
                
                for j, kx in enumerate(kx_um):
                    for k, ky in enumerate(ky_um):
                        results = self._calculate_single_point(E, kx, ky, params)
                        
                        Rs[global_i, j, k] = results['Rs']
                        Rp[global_i, j, k] = results['Rp']
                        S0[global_i, j, k] = results['S0']
                        S1[global_i, j, k] = results['S1']
                        S2[global_i, j, k] = results['S2']
                        S3[global_i, j, k] = results['S3']
                        phi[global_i, j, k] = results['phi']
                        chi[global_i, j, k] = results['chi']
        
        computation_time = time.time() - start_time
        memory_usage = self.estimate_memory_usage(params)
        
        return JonesResults(
            energy_eV=energy_eV,
            kx_um=kx_um,
            ky_um=ky_um,
            Rs=Rs, Rp=Rp,
            S0=S0, S1=S1, S2=S2, S3=S3,
            phi=phi, chi=chi,
            calculation_type="3D E-kx-ky mapping (chunked)",
            grid_resolution=params.resolution.value,
            computation_time=computation_time,
            memory_usage_mb=memory_usage
        )
    
    def _calculate_single_point(self, energy_eV: float, kx_um: float, ky_um: float,
                              params: CalculationParams) -> Dict[str, float]:
        """Calculate Jones matrix results for a single (E, kx, ky) point"""
        # Convert units
        w_spec = energy_eV * 8065.54429  # eV to cm⁻¹
        omega = 2 * np.pi * Constants.C * w_spec * 100  # rad/s
        kx = kx_um * 1e6  # μm⁻¹ to m⁻¹
        ky = ky_um * 1e6  # μm⁻¹ to m⁻¹
        k_par = np.hypot(kx, ky)
        
        # Build Jones matrix for multilayer structure
        G = self._build_jones_matrix(energy_eV, kx, ky, omega, k_par)
        
        # Input polarization
        delta = np.deg2rad(params.phase_diff_deg)
        Ein_lab = np.array([-params.ep_amplitude * np.exp(1j * delta), 
                           params.es_amplitude, 0])
        
        # Transform to local coordinates
        R = PhysicsUtils.local_basis_rotation(kx, ky)
        Ein_loc = R.T @ Ein_lab
        
        # Apply Jones matrix
        if G is not None:
            Eout_loc = G @ Ein_loc[:2]
            Eref_lab = R @ np.array([Eout_loc[0], Eout_loc[1], 0])
            Ex, Ey = Eref_lab[0], Eref_lab[1]
            
            # Calculate optical quantities
            S0, S1, S2, S3 = PhysicsUtils.stokes_parameters(Ex, Ey)
            phi, chi = PhysicsUtils.polarization_ellipse_parameters(S0, S1, S2, S3)
            
            Rs = abs(G[0, 0])**2
            Rp = abs(G[1, 1])**2
        else:
            # Fallback values
            Rs = Rp = S0 = S1 = S2 = S3 = phi = chi = 0.0
        
        return {
            'Rs': Rs, 'Rp': Rp,
            'S0': S0, 'S1': S1, 'S2': S2, 'S3': S3,
            'phi': phi, 'chi': chi
        }
    
    def _build_jones_matrix(self, energy_eV: float, kx: float, ky: float,
                          omega: float, k_par: float) -> Optional[np.ndarray]:
        """Build Jones matrix for multilayer structure"""
        # This is a simplified version - we'll expand this to use the full TMM formalism
        # For now, return identity matrix as placeholder
        return np.eye(2, dtype=complex)
    
    def calculate_kspace_slice(self, energy_eV: float, 
                             params: Optional[CalculationParams] = None) -> JonesResults:
        """
        Calculate 2D kx-ky cross-section at fixed energy.
        
        Parameters
        ----------
        energy_eV : float
            Fixed energy in eV
        params : CalculationParams, optional
            Calculation parameters
            
        Returns
        -------
        JonesResults
            2D k-space cross-section results
        """
        if params is None:
            params = self.default_params
        
        # Create single-energy parameters
        slice_params = CalculationParams(
            energy_start_eV=energy_eV,
            energy_stop_eV=energy_eV,
            energy_points=1,
            kx_start_um=params.kx_start_um,
            kx_stop_um=params.kx_stop_um,
            kx_points=params.kx_points,
            ky_start_um=params.ky_start_um,
            ky_stop_um=params.ky_stop_um,
            ky_points=params.ky_points,
            es_amplitude=params.es_amplitude,
            ep_amplitude=params.ep_amplitude,
            phase_diff_deg=params.phase_diff_deg,
            resolution=params.resolution
        )
        
        # Calculate 3D with single energy point
        results_3d = self.calculate_3d_mapping(slice_params)
        
        # Extract 2D slice
        return JonesResults(
            energy_eV=results_3d.energy_eV,
            kx_um=results_3d.kx_um,
            ky_um=results_3d.ky_um,
            Rs=results_3d.Rs[0, :, :],
            Rp=results_3d.Rp[0, :, :],
            S0=results_3d.S0[0, :, :],
            S1=results_3d.S1[0, :, :],
            S2=results_3d.S2[0, :, :],
            S3=results_3d.S3[0, :, :],
            phi=results_3d.phi[0, :, :],
            chi=results_3d.chi[0, :, :],
            calculation_type=f"2D kx-ky slice at {energy_eV:.2f} eV",
            grid_resolution=results_3d.grid_resolution,
            computation_time=results_3d.computation_time,
            memory_usage_mb=results_3d.memory_usage_mb
        )
    
    def calculate_dispersion_line(self, k_direction: str, k_value: float,
                                params: Optional[CalculationParams] = None) -> JonesResults:
        """
        Calculate 1D dispersion along k-direction.
        
        Parameters
        ----------
        k_direction : str
            Direction to fix ('kx' or 'ky')
        k_value : float
            Fixed k-value in μm⁻¹
        params : CalculationParams, optional
            Calculation parameters
            
        Returns
        -------
        JonesResults
            1D dispersion results
        """
        if params is None:
            params = self.default_params
        
        # Create dispersion parameters
        if k_direction == 'kx':
            disp_params = CalculationParams(
                energy_start_eV=params.energy_start_eV,
                energy_stop_eV=params.energy_stop_eV,
                energy_points=params.energy_points,
                kx_start_um=k_value,
                kx_stop_um=k_value,
                kx_points=1,
                ky_start_um=params.ky_start_um,
                ky_stop_um=params.ky_stop_um,
                ky_points=params.ky_points,
                es_amplitude=params.es_amplitude,
                ep_amplitude=params.ep_amplitude,
                phase_diff_deg=params.phase_diff_deg,
                resolution=params.resolution
            )
        else:  # k_direction == 'ky'
            disp_params = CalculationParams(
                energy_start_eV=params.energy_start_eV,
                energy_stop_eV=params.energy_stop_eV,
                energy_points=params.energy_points,
                kx_start_um=params.kx_start_um,
                kx_stop_um=params.kx_stop_um,
                kx_points=params.kx_points,
                ky_start_um=k_value,
                ky_stop_um=k_value,
                ky_points=1,
                es_amplitude=params.es_amplitude,
                ep_amplitude=params.ep_amplitude,
                phase_diff_deg=params.phase_diff_deg,
                resolution=params.resolution
            )
        
        # Calculate 3D with fixed k-direction
        results_3d = self.calculate_3d_mapping(disp_params)
        
        # Extract 1D dispersion
        if k_direction == 'kx':
            Rs = results_3d.Rs[:, 0, :]
            Rp = results_3d.Rp[:, 0, :]
            S0 = results_3d.S0[:, 0, :]
            S1 = results_3d.S1[:, 0, :]
            S2 = results_3d.S2[:, 0, :]
            S3 = results_3d.S3[:, 0, :]
            phi = results_3d.phi[:, 0, :]
            chi = results_3d.chi[:, 0, :]
        else:
            Rs = results_3d.Rs[:, :, 0]
            Rp = results_3d.Rp[:, :, 0]
            S0 = results_3d.S0[:, :, 0]
            S1 = results_3d.S1[:, :, 0]
            S2 = results_3d.S2[:, :, 0]
            S3 = results_3d.S3[:, :, 0]
            phi = results_3d.phi[:, :, 0]
            chi = results_3d.chi[:, :, 0]
        
        return JonesResults(
            energy_eV=results_3d.energy_eV,
            kx_um=results_3d.kx_um,
            ky_um=results_3d.ky_um,
            Rs=Rs, Rp=Rp,
            S0=S0, S1=S1, S2=S2, S3=S3,
            phi=phi, chi=chi,
            calculation_type=f"1D dispersion with {k_direction}={k_value:.2f} μm⁻¹",
            grid_resolution=results_3d.grid_resolution,
            computation_time=results_3d.computation_time,
            memory_usage_mb=results_3d.memory_usage_mb
        ) 