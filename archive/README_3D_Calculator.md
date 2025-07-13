# 3D TMM Calculator - E-kx-ky Reflectance Analysis

## Overview

The 3D TMM (Transfer Matrix Method) Calculator extends the functionality of the original 2D codes (`Jones_disperion.py` and `Jones_Xect.py`) to calculate reflectance and Stokes parameters in full 3D space: **Energy-kx-ky**.

### Original 2D Codes:
- `Jones_disperion.py`: Calculates E-k dispersion (Energy vs single k-vector component)
- `Jones_Xect.py`: Calculates kx-ky cross-section (kx vs ky at fixed energy)

### New 3D Code:
- `Jones_3D_Calculator.py`: Calculates E-kx-ky (Energy vs kx vs ky simultaneously)

## Features

### Core Functionality
- **3D Calculation**: Simultaneous calculation across Energy-kx-ky space
- **Jones Matrix Method**: Full polarization analysis using Jones calculus
- **Stokes Parameters**: Complete polarimetric analysis (S0, S1, S2, S3, φ, χ)
- **Reflectance Components**: Separate s-polarized (Rs) and p-polarized (Rp) reflectance
- **Layered Structures**: Support for arbitrary multilayer structures
- **Material Database**: Built-in materials (Air, SiO2) with custom material file loading

### GUI Features
- **Interactive Interface**: User-friendly GUI for parameter setting
- **Layer Management**: Add/edit/delete layers with thickness and material selection
- **Grid Presets**: Coarse (testing), Medium, and Fine (production) grid options
- **Progress Monitoring**: Real-time progress bar and status updates
- **Results Visualization**: Built-in 3D visualization tools
- **Data Export**: Save results in numpy format for further analysis

### Performance Considerations
- **Coarse Grid**: ~10-100 points for quick testing
- **Medium Grid**: ~1,000-10,000 points for intermediate analysis
- **Fine Grid**: ~100,000+ points for production results
- **Computational Warning**: Automatically warns for large calculations

## Installation & Dependencies

### Required Libraries
```bash
pip install numpy matplotlib tkinter
```

### File Structure
```
TMM/archive/
├── Jones_3D_Calculator.py     # Main 3D calculator
├── test_3d_calculator.py      # Test script
├── Jones_disperion.py         # Original E-k dispersion code
├── Jones_Xect.py              # Original kx-ky cross-section code
└── README_3D_Calculator.md    # This file
```

## Usage

### 1. GUI Application
```python
python Jones_3D_Calculator.py
```

**Step-by-step GUI usage:**
1. **Set Energy Range**: Enter E_start, E_stop, and dE (step size)
2. **Set kx Range**: Enter kx_start, kx_stop, and dkx
3. **Set ky Range**: Enter ky_start, ky_stop, and dky
4. **Choose Grid Preset**: 
   - "Coarse (Test)": For quick testing
   - "Medium": For intermediate analysis
   - "Fine (Final)": For production results
5. **Set Input Polarization**: |Es|, |Ep|, and phase difference Δφ
6. **Configure Layers**: Add/edit/delete layers using the layer table
7. **Load Materials**: Use "Load Material File" for custom materials
8. **Calculate**: Click "Calculate 3D" to start computation
9. **Visualize**: Use "Visualize Results" to view 3D plots
10. **Save**: Export results using "Save Results"

### 2. Command Line Usage
```python
from Jones_3D_Calculator import calc_3d_reflectance

# Define layer structure
layers = [
    ('Air', 1e6),    # Semi-infinite substrate
    ('SiO2', 500),   # 500 nm SiO2 layer
    ('Air', 1e6)     # Semi-infinite superstrate
]

# Calculate 3D reflectance
results = calc_3d_reflectance(
    layers,                    # Layer structure
    E_start=1.0, E_stop=3.0, dE=0.1,    # Energy range
    kx_start=-10, kx_stop=10, dkx=0.5,   # kx range
    ky_start=-10, ky_stop=10, dky=0.5,   # ky range
    a=1.0, b=0.0, delta_phi_deg=0.0      # Input polarization
)

E_vals, kx_vals, ky_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = results
```

### 3. Test Script
```python
python test_3d_calculator.py
```

Choose between:
1. Command-line calculation test (demonstrates basic usage)
2. GUI application test (launches full interface)

## Output Parameters

### Reflectance Components
- **Rs**: s-polarized reflectance (|r_s|²)
- **Rp**: p-polarized reflectance (|r_p|²)

### Stokes Parameters
- **S0**: Total intensity (|Ex|² + |Ey|²)
- **S1**: Linear polarization difference (|Ex|² - |Ey|²)/S0
- **S2**: Linear polarization at 45° (2·Re(Ex·Ey*)/S0)
- **S3**: Circular polarization (2·Im(Ex·Ey*)/S0)

### Polarization Angles
- **φ (phi)**: Polarization angle (0.5·arctan2(S2, S1))
- **χ (chi)**: Ellipticity angle (0.5·arcsin(S3/S0))

## Layer Structure Format

### Built-in Materials
- **'Air'**: Vacuum/air (n=1.0)
- **'SiO2'**: Silicon dioxide (n=1.5)

### Custom Materials
Material data files should contain columns:
```
frequency(cm⁻¹)  Re(εx)  Im(εx)  Re(εy)  Im(εy)  Re(εz)  Im(εz)
```

### Layer Definition
```python
layers = [
    ('Material1', thickness1_nm),
    ('Material2', thickness2_nm),
    ...
]
```

## Coordinate System

### Laboratory Frame
- **x-axis**: In-plane direction
- **y-axis**: In-plane direction  
- **z-axis**: Out-of-plane (growth direction)

### Wave Vector
- **kx**: x-component of wave vector (μm⁻¹)
- **ky**: y-component of wave vector (μm⁻¹)
- **kz**: z-component (calculated from dispersion)

### Polarization Convention
- **s-polarization**: Electric field perpendicular to incident plane
- **p-polarization**: Electric field parallel to incident plane

## Performance Guidelines

### Grid Size Recommendations

| Application | Energy Points | kx Points | ky Points | Total Points | Est. Time |
|-------------|---------------|-----------|-----------|--------------|-----------|
| Quick Test  | 5-10         | 10-20     | 10-20     | 500-4,000   | Seconds   |
| Development | 20-50        | 20-50     | 20-50     | 8,000-125,000| Minutes   |
| Production  | 50-200       | 50-200    | 50-200    | 125,000-8M  | Hours     |

### Memory Requirements
- **Coarse**: ~10 MB RAM
- **Medium**: ~100 MB RAM  
- **Fine**: ~1-10 GB RAM

### Computational Complexity
- **Scaling**: O(nE × nkx × nky × nlayers)
- **Bottleneck**: Jones matrix calculations for each (E, kx, ky) point
- **Parallelization**: Currently single-threaded (future enhancement)

## Output File Formats

### Coordinate Arrays
- `E_vals.dat`: Energy values (eV)
- `kx_vals.dat`: kx values (μm⁻¹)
- `ky_vals.dat`: ky values (μm⁻¹)

### 3D Data Arrays (Binary)
- `Rs.npy`: s-polarized reflectance [nE × nkx × nky]
- `Rp.npy`: p-polarized reflectance [nE × nkx × nky]
- `S0.npy`: Total intensity [nE × nkx × nky]
- `S1.npy`: Linear polarization S1 [nE × nkx × nky]
- `S2.npy`: Linear polarization S2 [nE × nkx × nky]
- `S3.npy`: Circular polarization [nE × nkx × nky]
- `phi.npy`: Polarization angle [nE × nkx × nky]
- `chi.npy`: Ellipticity angle [nE × nkx × nky]

### Metadata
- `metadata.txt`: Calculation parameters and layer information

## Comparison with Original 2D Codes

| Feature | Jones_disperion.py | Jones_Xect.py | Jones_3D_Calculator.py |
|---------|-------------------|---------------|------------------------|
| **Calculation Space** | E-k (2D) | kx-ky (2D) | E-kx-ky (3D) |
| **Fixed Parameter** | kx or ky | Energy | None |
| **Computational Time** | Minutes | Minutes | Hours |
| **Memory Usage** | Low | Low | High |
| **Visualization** | 2D plots | 2D plots | 3D visualization |
| **Data Export** | Text files | Text files | Binary + text |

## Applications

### Research Applications
- **Metamaterial Design**: Analyze dispersion properties
- **Photonic Crystals**: Study band structure and polarization
- **Thin Film Optics**: Optimize multilayer coatings
- **Plasmonics**: Investigate surface plasmon modes
- **Quantum Optics**: Analyze cavity QED systems

### Industrial Applications
- **Optical Coating Design**: AR/HR coatings optimization
- **Solar Cell Optimization**: Light trapping structures
- **Display Technology**: Polarization management
- **Sensor Development**: Optical sensing applications

## Troubleshooting

### Common Issues
1. **Memory Error**: Reduce grid size or use coarse grid
2. **Calculation Stuck**: Check for proper layer structure
3. **Visualization Slow**: Reduce data size for plotting
4. **Material Loading**: Verify file format and path

### Performance Tips
1. **Start with Coarse Grid**: Always test with small grids first
2. **Monitor Progress**: Use progress callback for long calculations
3. **Save Regularly**: Export intermediate results
4. **Check Layer Structure**: Ensure proper material definitions

## Future Enhancements

### Planned Features
- **Parallel Processing**: Multi-core CPU utilization
- **GPU Acceleration**: CUDA/OpenCL support
- **Advanced Visualization**: Interactive 3D plots
- **Batch Processing**: Multiple structure calculations
- **Material Database**: Extended material library

### Optimization Opportunities
- **Memory Management**: Reduce RAM usage for large grids
- **Numerical Stability**: Improve convergence for extreme parameters
- **User Interface**: Enhanced GUI with real-time previews

## References

1. Transfer Matrix Method in optics
2. Jones calculus for polarization analysis
3. Stokes parameters in polarimetry
4. Multilayer thin film optics

## License

This code is provided for academic and research purposes. Please cite appropriately when used in publications.

---

**Contact**: For questions or suggestions, please refer to the original codebase documentation or contact the development team. 