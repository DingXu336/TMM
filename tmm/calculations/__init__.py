"""
TMM calculations module containing specialized calculators.

This module provides high-level calculation interfaces:
- TMMCalculator: General TMM calculations
- ReflectivityCalculator: Reflectivity and transmission calculations
- DispersionCalculator: Energy/wavelength dispersion calculations
- PolarizationAnalyzer: Jones matrix and Stokes parameter analysis
"""

from .tmm_calculator import TMMCalculator
from .reflectivity_calculator import ReflectivityCalculator
from .dispersion_calculator import DispersionCalculator
from .polarization_analyzer import PolarizationAnalyzer

__all__ = [
    "TMMCalculator",
    "ReflectivityCalculator", 
    "DispersionCalculator",
    "PolarizationAnalyzer",
] 