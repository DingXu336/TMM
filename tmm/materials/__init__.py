"""
Materials module for material property management.

This module provides:
- Material: Individual material representation
- MaterialLoader: Loading material data from files
- Built-in materials: Common materials like Air, SiO2, etc.
"""

from .material import Material
from .material_loader import MaterialLoader
from .builtin_materials import get_builtin_material

__all__ = ["Material", "MaterialLoader", "get_builtin_material"] 