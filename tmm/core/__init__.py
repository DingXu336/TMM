"""
Core TMM module containing fundamental classes and data structures.

This module provides the basic building blocks for TMM calculations:
- TMM: Main Transfer Matrix Method class
- Layer: Individual layer representation
- MaterialDatabase: Material property management
"""

# Import only implemented modules
try:
    from .layer import Layer, create_substrate, create_superstrate, create_film
    __all__ = ["Layer", "create_substrate", "create_superstrate", "create_film"]
except ImportError:
    __all__ = [] 