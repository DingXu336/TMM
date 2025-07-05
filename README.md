# TMM (Transfer Matrix Method) Package

A comprehensive Python package for Transfer Matrix Method calculations in optical layered structures, designed for physical chemistry research.

**Author**: Ding Xu (Physical Chemistry Researcher)  
**Version**: 1.0.0  
**License**: MIT  

## Overview

The TMM package provides a complete suite of tools for calculating optical properties of layered structures using the Transfer Matrix Method. It's designed with modern software engineering principles and provides both programmatic APIs and graphical user interfaces for interactive use.

### Key Features

- **Comprehensive TMM calculations**: Reflectivity, transmission, and absorption
- **Jones matrix formalism**: Full polarization analysis with Stokes parameters
- **Material database**: Built-in materials and custom material loading
- **2D k-space calculations**: Angle-resolved optical properties
- **Dispersion analysis**: Energy and wavelength-dependent calculations
- **Professional APIs**: Well-documented, type-hinted interfaces
- **GUI applications**: Interactive interfaces for calculations
- **Extensive validation**: Input validation and error handling
- **Scientific accuracy**: Based on established electromagnetic theory

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/dingxu/tmm-optics.git
cd tmm-optics

# Install in development mode
pip install -e .

# Install with GUI dependencies
pip install -e ".[gui]"

# Install with development dependencies
pip install -e ".[dev]"
```

### Requirements

- Python ≥ 3.8
- NumPy ≥ 1.20.0
- SciPy ≥ 1.7.0
- Matplotlib ≥ 3.5.0
- Optional: tkinter (for GUI, usually included with Python)

## Quick Start

### Basic TMM Calculation

```python
import numpy as np
from tmm import TMM, Layer, TMMCalculator
from tmm.materials import get_builtin_material

# Create materials
air = get_builtin_material("Air", energy_eV=2.0)
sio2 = get_builtin_material("SiO2", energy_eV=2.0)

# Create layer structure
layers = [
    Layer(air, thickness=np.inf),      # Semi-infinite superstrate
    Layer(sio2, thickness=500e-9),     # 500 nm SiO2 film
    Layer(air, thickness=np.inf)       # Semi-infinite substrate
]

# Create TMM calculator
calculator = TMMCalculator(layers)

# Calculate reflectivity
energy_eV = np.linspace(1.0, 3.0, 100)
results = calculator.calculate_reflectivity(
    energy_eV=energy_eV,
    angle_deg=0.0,
    polarization="s"
)

print(f"Reflectivity at 2.0 eV: {results['reflectivity'][50]:.3f}")
```

### Material Loading

```python
from tmm.materials import MaterialLoader

# Load material from n,k data file
material = MaterialLoader.load_from_nk_file(
    "materials/Au_JC.txt",
    material_name="Gold",
    wavelength_unit="nm"
)

# Load from dielectric function file
material = MaterialLoader.load_from_dielectric_file(
    "materials/Si_dielectric.txt",
    energy_unit="eV"
)
```

### 2D k-space Calculations

```python
from tmm.calculations import ReflectivityCalculator

# Create calculator
calculator = ReflectivityCalculator(layers)

# Calculate 2D reflectivity map
kx_range = np.linspace(0, 25, 100)  # μm⁻¹
ky_range = np.linspace(0, 25, 100)  # μm⁻¹

results = calculator.calculate_2d_reflectivity(
    energy_eV=2.0,
    kx_range=kx_range,
    ky_range=ky_range,
    polarization_s_amplitude=1.0,
    polarization_p_amplitude=0.0
)

# Results contain reflectivity, Stokes parameters, etc.
```

### Dispersion Analysis

```python
from tmm.calculations import DispersionCalculator

# Create calculator
calculator = DispersionCalculator(layers)

# Calculate energy dispersion
energy_range = np.linspace(1.0, 4.0, 100)
k_parallel = 10.0  # μm⁻¹

results = calculator.calculate_energy_dispersion(
    energy_range=energy_range,
    k_parallel=k_parallel,
    polarization="s"
)
```

## GUI Applications

The package includes several GUI applications for interactive calculations:

```bash
# Launch main TMM GUI
tmm-gui

# Convert material files
tmm-convert
```

### GUI Features

- **Layer structure editor**: Visual layer structure management
- **Material database browser**: Browse and load materials
- **Real-time calculations**: Interactive parameter adjustment
- **2D visualization**: Angle-resolved reflectivity maps
- **Dispersion plots**: Energy and wavelength dependence
- **Export capabilities**: Save results and plots

## Scientific Background

### Transfer Matrix Method

The Transfer Matrix Method (TMM) is a computational technique for calculating optical properties of layered structures. It's based on the continuity of tangential electric and magnetic fields at interfaces between layers.

**Key Equations:**

For s-polarized light (TE):
```
kz = √(ε * k₀² - k‖²)
```

For p-polarized light (TM):
```
kz = √(εxx * k₀² - (εxx/εzz) * k‖²)
```

### Jones Matrix Formalism

The package implements full Jones matrix calculations for polarization analysis:

```
J = (Y₂ - Y₁) / (Y₂ + Y₁)
```

Where Y₁ and Y₂ are the optical admittances of adjacent layers.

### Stokes Parameters

Complete polarization state analysis using Stokes parameters:

- S₀: Total intensity
- S₁: Linear polarization (horizontal - vertical)
- S₂: Linear polarization (45° - 135°)
- S₃: Circular polarization

## API Reference

### Core Classes

- **`TMM`**: Main TMM calculation engine
- **`Layer`**: Individual layer representation
- **`Material`**: Material property management
- **`MaterialDatabase`**: Material database management

### Calculators

- **`TMMCalculator`**: General TMM calculations
- **`ReflectivityCalculator`**: Reflectivity and transmission
- **`DispersionCalculator`**: Energy/wavelength dispersion
- **`PolarizationAnalyzer`**: Jones matrix and Stokes analysis

### Utilities

- **`Constants`**: Physical constants
- **`PhysicsUtils`**: Physics calculations and transformations
- **`ValidationUtils`**: Input validation and error checking

## Examples and Tutorials

The `examples/` directory contains comprehensive examples:

- **`basic_tmm.py`**: Basic TMM calculations
- **`material_loading.py`**: Loading custom materials
- **`2d_reflectivity.py`**: 2D k-space calculations
- **`dispersion_analysis.py`**: Dispersion calculations
- **`polarization_analysis.py`**: Jones matrix calculations
- **`fabry_perot.py`**: Fabry-Pérot cavity analysis
- **`multilayer_coating.py`**: Antireflection coating design

### Jupyter Notebooks

Interactive tutorials available in `tutorials/`:

- **`01_introduction.ipynb`**: Package overview
- **`02_materials.ipynb`**: Material system
- **`03_basic_calculations.ipynb`**: Basic TMM calculations
- **`04_advanced_features.ipynb`**: Advanced features
- **`05_polarization.ipynb`**: Polarization analysis

## Testing

Run the test suite:

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=tmm

# Run specific test category
pytest tests/test_core.py
pytest tests/test_materials.py
pytest tests/test_calculations.py
```

## Performance Considerations

- **Vectorization**: All calculations are vectorized using NumPy
- **Memory management**: Efficient memory usage for large calculations
- **Parallel computing**: Support for parallel calculations (future)
- **Optimization**: Optimized algorithms for common use cases

## Contributing

Contributions are welcome! Please see `CONTRIBUTING.md` for guidelines.

### Development Setup

```bash
# Clone repository
git clone https://github.com/dingxu/tmm-optics.git
cd tmm-optics

# Install development dependencies
pip install -e ".[dev]"

# Run code formatting
black tmm/
flake8 tmm/

# Run type checking
mypy tmm/
```

## Citation

If you use this package in your research, please cite:

```bibtex
@software{tmm_optics,
  author = {Xu, Ding},
  title = {TMM: Transfer Matrix Method for Optical Calculations},
  url = {https://github.com/dingxu/tmm-optics},
  version = {1.0.0},
  year = {2024}
}
```

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

## Acknowledgments

- Based on electromagnetic theory and TMM formalism
- Inspired by research in physical chemistry and optics
- Thanks to the scientific Python community

## Support

- **Documentation**: [https://tmm-optics.readthedocs.io](https://tmm-optics.readthedocs.io)
- **Issues**: [https://github.com/dingxu/tmm-optics/issues](https://github.com/dingxu/tmm-optics/issues)
- **Discussions**: [https://github.com/dingxu/tmm-optics/discussions](https://github.com/dingxu/tmm-optics/discussions)

## Archive

The `archive/` directory contains the original Python files from before refactoring. These files are deprecated and preserved for historical reference only. **Use the new `tmm` package for all new work.**

See `archive/README.md` for details about the migration from original files to the new structure.

## Changelog

### Version 1.0.0 (2024)

- Initial release
- Complete TMM implementation
- Material database system
- GUI applications
- Comprehensive documentation
- Full test coverage
- Refactored from original standalone scripts

---

For more information, please refer to the [documentation](https://tmm-optics.readthedocs.io) or contact the author.
