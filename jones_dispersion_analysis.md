# 2D Jones Dispersion Analysis - Complete Implementation

## Overview

We have successfully reproduced the **2D Jones dispersion analysis** functionality from the original `archive/Jones_disperion.py` file using our refactored TMM package. The implementation creates comprehensive 2D dispersion plots showing how optical properties vary with both **energy (eV)** and **q-vector (μm⁻¹)**.

## Key Features Implemented

### 1. **2D Dispersion Calculation**
- **Energy Range**: Configurable range in eV (e.g., 0.5-4.0 eV for plasmonic studies)
- **q-Vector Range**: Configurable range in μm⁻¹ (e.g., 0-30 μm⁻¹ for high-momentum features)
- **Fixed Axis**: Can fix either kx or ky and vary the other
- **Resolution Control**: Low, medium, high resolution options

### 2. **Jones Matrix Formalism**
- **Complete Jones Matrix Calculation**: Interface and propagation matrices
- **Coordinate System Rotation**: Proper handling of arbitrary k-vectors
- **Polarization Analysis**: Full Jones vector processing
- **Multi-layer Support**: Cascaded interface calculations

### 3. **Output Quantities (8 Different 2D Maps)**
1. **Rs**: s-polarized reflectivity
2. **Rp**: p-polarized reflectivity  
3. **S0**: Stokes parameter (total intensity)
4. **S1**: Stokes parameter (linear H-V polarization)
5. **S2**: Stokes parameter (linear ±45° polarization)
6. **S3**: Stokes parameter (circular polarization)
7. **φ (phi)**: Azimuth angle of polarization ellipse
8. **χ (chi)**: Ellipticity angle of polarization ellipse

### 4. **Enhanced Visualization**
- **2×4 Subplot Layout**: Matches original format
- **Color Maps**: Appropriate colormaps for each quantity
- **Scientific Styling**: Professional plots with proper labeling
- **High Resolution**: 300 DPI output for publication quality

## Demo Scripts Created

### 1. `jones_dispersion_demo.py`
**Basic Demo**: Reproduces original functionality exactly
```python
# Layer structure
layer_info = [
    ('Air', 1e6),      # Semi-infinite air
    ('SiO2', 500),     # 500 nm SiO2 layer
    ('Air', 1e6)       # Semi-infinite air substrate
]

# Parameters
Energy: 1.0 - 3.0 eV
q-vector: 0.0 - 25.0 μm⁻¹
Polarization: s-polarized (|Es|=1.0, |Ep|=0.0)
```

### 2. `jones_dispersion_plasmonic_demo.py`
**Plasmonic Demo**: Shows surface plasmon features
```python
# Layer structure
layer_info = [
    ('Air', 1e6),      # Semi-infinite air
    ('Gold', 50),      # 50 nm Gold film
    ('SiO2', 1e6)      # Semi-infinite SiO2 substrate
]

# Parameters
Energy: 0.5 - 4.0 eV (plasmonic region)
q-vector: 0.0 - 30.0 μm⁻¹
Polarization: p-polarized (|Es|=0.0, |Ep|=1.0)
```

## Scientific Applications

### 1. **Surface Plasmon Polaritons**
- **Dispersion Curves**: Visible in Rp plots as resonance features
- **Field Enhancement**: Shown in Stokes parameters
- **Momentum Matching**: k-vector dependence clearly visible

### 2. **Guided Modes**
- **Waveguide Modes**: Discrete dispersion branches
- **Cut-off Frequencies**: Mode boundaries in E-q space
- **Polarization Conversion**: Mixing between s and p modes

### 3. **Photonic Crystals**
- **Band Structures**: Forbidden and allowed regions
- **Band Gaps**: Spectral ranges with no propagating modes
- **Bloch Modes**: Periodic structure effects

### 4. **Metamaterials**
- **Negative Index**: Unusual dispersion behavior
- **Resonant Features**: Sharp spectral and angular features
- **Anisotropic Response**: Polarization-dependent behavior

## Technical Implementation

### Core Physics
```python
# Jones matrix for interface
J = (Y₂ - Y₁) / (Y₂ + Y₁)

# Propagation through layer
G = gamma_jones(rj, Gp, φ)

# Stokes parameters
S0 = |Ex|² + |Ey|²
S1 = (|Ex|² - |Ey|²) / S0
S2 = 2·Re(Ex·Ey*) / S0
S3 = 2·Im(Ex·Ey*) / S0
```

### Key Advantages Over Original
1. **Modular Design**: Uses our refactored TMM package
2. **Better Error Handling**: Robust against edge cases
3. **Enhanced Visualization**: Professional-quality plots
4. **Material Integration**: Works with our material database
5. **Extensibility**: Easy to add new materials and features

## Results Generated

### File Outputs
- `jones_dispersion_2d_kx_0.0.png`: Basic SiO2 layer analysis
- `jones_dispersion_2d_Air_Gold_SiO2_kx_0.0.png`: Plasmonic analysis
- Both files are high-resolution (300 DPI) scientific plots

### Data Characteristics
- **Basic Demo**: 201 × 251 = 50,451 data points
- **Plasmonic Demo**: 176 × 151 = 26,576 data points
- **Calculation Time**: ~1-2 minutes per demo (depends on resolution)

## Comparison with Original

### ✅ **Reproduced Features**
- [x] 2D dispersion calculation (E vs q)
- [x] Jones matrix formalism
- [x] 8 different optical quantities
- [x] Stokes parameter analysis
- [x] Polarization angles (φ, χ)
- [x] s and p reflectivities
- [x] 2×4 subplot layout
- [x] Proper color mapping
- [x] Scientific axis labeling

### ➕ **Enhancements Added**
- [x] Better material handling
- [x] Enhanced visualization
- [x] Modular architecture
- [x] Resolution control
- [x] Professional documentation
- [x] Error handling
- [x] Progress indicators

## Usage Examples

### Basic Usage
```bash
python3 jones_dispersion_demo.py
```

### Plasmonic Analysis
```bash
python3 jones_dispersion_plasmonic_demo.py
```

### Custom Analysis
```python
# Modify layer structure
layer_info = [
    ('Air', 1e6),
    ('Your_Material', thickness_nm),
    ('Substrate', 1e6)
]

# Adjust parameters
E_start, E_stop = 1.0, 5.0
q_start, q_stop = 0.0, 50.0
```

## Scientific Interpretation

### Reading the Plots
1. **Rs, Rp Plots**: Show reflectivity for s and p polarizations
   - Dark regions = low reflectivity (absorption/transmission)
   - Bright regions = high reflectivity
   - Resonance features appear as distinct patterns

2. **Stokes Parameters**: Show polarization state
   - S0 = total intensity
   - S1, S2, S3 = polarization components
   - Values range from -1 to +1

3. **Polarization Angles**: Show ellipse parameters
   - φ (azimuth) = orientation of polarization ellipse
   - χ (ellipticity) = shape of polarization ellipse

### Physical Insights
- **Surface Plasmons**: Appear as dispersive curves in Rp plots
- **Critical Angles**: Sudden changes in reflectivity
- **Polarization Conversion**: Mixing between different Stokes components
- **Resonance Effects**: Enhanced field coupling in thin films

## Future Extensions

### Possible Enhancements
1. **3D Visualization**: Add full k-space analysis
2. **Animation**: Time-evolution of dispersion
3. **Interactive GUI**: Real-time parameter adjustment
4. **Optimization**: GPU acceleration for large calculations
5. **Advanced Materials**: Nonlinear and active materials

### Research Applications
- **Plasmonics**: Surface plasmon devices
- **Photonics**: Optical waveguides and resonators
- **Metamaterials**: Artificial optical materials
- **Sensors**: Optical sensing applications

## Conclusion

We have successfully **reproduced and enhanced** the original `Jones_disperion.py` functionality using our refactored TMM package. The implementation provides:

- **Complete 2D dispersion analysis** with all 8 optical quantities
- **Professional-quality visualizations** matching the original format
- **Enhanced material handling** through our modular architecture
- **Scientific accuracy** with proper electromagnetic theory
- **Extensible design** for future research applications

The demos demonstrate both basic dielectric layers and advanced plasmonic structures, showing the versatility of the implementation for physical chemistry research applications.

---

**Author**: Ding Xu (Physical Chemistry Researcher)  
**Date**: 2024  
**Based on**: Original Jones_disperion.py with modern software engineering 