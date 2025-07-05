"""
GUI module for interactive TMM applications.

This module provides graphical user interfaces for:
- 2D reflectivity calculations
- Dispersion analysis
- Material property visualization
- Jones matrix analysis
"""

from .main import main
from .tmm_gui import TMMGui
from .dispersion_gui import DispersionGui

__all__ = ["main", "TMMGui", "DispersionGui"] 