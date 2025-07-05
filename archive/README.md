# Archive - Original TMM Files

This directory contains the original Python files from the TMM project before refactoring.

## Archived Files

- **`TMM_2D.py`** - Original main TMM calculation script with integrated GUI
- **`Jones_Xect.py`** - Original Jones matrix cross-section calculations
- **`Jones_disperion.py`** - Original Jones matrix dispersion analysis 
- **`Jones_dispersion_unfinished.py`** - Incomplete dispersion implementation
- **`eps_convert_from_ref-web.py`** - Original material property conversion utility

## Migration Status

These files have been **completely refactored** into the new modular package structure.

### Original → New Structure Mapping

| Original File | New Location |
|---------------|--------------|
| `TMM_2D.py` | → `tmm/core/tmm.py` + `tmm/gui/tmm_gui.py` |
| `Jones_Xect.py` | → `tmm/calculations/polarization_analyzer.py` |
| `Jones_disperion.py` | → `tmm/calculations/dispersion_calculator.py` |
| `eps_convert_from_ref-web.py` | → `tmm/utils/converter.py` + `tmm/materials/material_loader.py` |

## Important Notes

⚠️ **These files are deprecated and should not be used for new work.**

✅ **Use the new `tmm` package instead** for better code organization, documentation, and features.

## Preservation Purpose

These files are preserved for historical reference and migration assistance only.

---

**Recommendation**: Use the new `tmm` package for all future work. 