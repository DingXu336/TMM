"""
TMM calculations module containing specialized calculators.

This module provides high-level calculation interfaces:
- JonesUnifiedCalculator: Unified 3D E-kx-ky analysis
- Jones3DAnalyzer: Specialized 3D Jones matrix analyzer
"""

# Import unified Jones calculator components (the ones that actually exist)
from .jones_unified_calculator import JonesUnifiedCalculator, CalculationParams, JonesResults, GridResolution
from .jones_3d_analyzer import Jones3DAnalyzer

# TODO: Import other calculators when they are implemented
# from .tmm_calculator import TMMCalculator
# from .reflectivity_calculator import ReflectivityCalculator
# from .dispersion_calculator import DispersionCalculator
# from .polarization_analyzer import PolarizationAnalyzer

__all__ = [
    "JonesUnifiedCalculator",
    "CalculationParams", 
    "JonesResults",
    "GridResolution",
    "Jones3DAnalyzer",
    # TODO: Add other calculators when implemented
    # "TMMCalculator",
    # "ReflectivityCalculator", 
    # "DispersionCalculator",
    # "PolarizationAnalyzer",
] 