"""
TMM (Transfer Matrix Method) Package

A comprehensive Python package for Transfer Matrix Method calculations
in optical layered structures, designed for physical chemistry research.

This package provides:
- Core TMM calculations for reflectivity and transmission
- Jones matrix formalism for polarization analysis
- Material database management
- GUI applications for interactive calculations
- Utilities for data conversion and analysis

Author: Ding Xu (Physical Chemistry Researcher)
License: MIT
"""

# Core components (import only what exists)
try:
    from .core.layer import Layer, create_substrate, create_superstrate, create_film
except ImportError:
    Layer = None
    
try:
    from .materials.material import Material
    from .materials.material_loader import MaterialLoader  
    from .materials.builtin_materials import get_builtin_material
except ImportError:
    Material = None
    MaterialLoader = None
    get_builtin_material = None
    
try:
    from .utils.constants import Constants
    from .utils.physics_utils import PhysicsUtils
except ImportError:
    Constants = None
    PhysicsUtils = None

__version__ = "1.0.0"
__author__ = "Ding Xu"
__email__ = "dingxu@example.com"
__license__ = "MIT"

__all__ = [
    "TMM",
    "Layer", 
    "MaterialDatabase",
    "TMMCalculator",
    "ReflectivityCalculator", 
    "DispersionCalculator",
    "PolarizationAnalyzer",
    "Material",
    "MaterialLoader",
    "Constants",
    "PhysicsUtils",
] 