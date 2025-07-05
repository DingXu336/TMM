"""
GUI module for interactive TMM applications.

This module provides graphical user interfaces for:
- Unified 3D E-kx-ky analysis
"""

# Import unified Jones GUI (the one that actually exists)
from .jones_unified_gui import JonesUnifiedGUI

# TODO: Import other GUI components when they are implemented
# from .main import main
# from .tmm_gui import TMMGui
# from .dispersion_gui import DispersionGui

__all__ = ["JonesUnifiedGUI"] 