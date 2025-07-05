#!/usr/bin/env python3
"""
Material file converter utility for TMM package.

This utility provides command-line conversion between different
material file formats used in optical calculations.

Author: Ding Xu
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from ..materials.material_loader import MaterialLoader

def convert_file(input_path: str, 
                output_path: Optional[str] = None,
                input_format: Optional[str] = None,
                output_format: str = "dielectric",
                wavelength_unit: str = "um",
                energy_unit: str = "eV") -> None:
    """
    Convert material file between formats.
    
    Parameters
    ----------
    input_path : str
        Path to input file
    output_path : str, optional
        Path to output file (auto-generated if None)
    input_format : str, optional
        Input format ('nk', 'dielectric', 'json', auto-detect if None)
    output_format : str, default='dielectric'
        Output format ('nk', 'dielectric', 'json')
    wavelength_unit : str, default='um'
        Unit for wavelength in n,k files
    energy_unit : str, default='eV'
        Unit for energy in dielectric files
    """
    input_path = Path(input_path)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    # Auto-detect input format if not specified
    if input_format is None:
        input_format = MaterialLoader.detect_file_format(input_path)
        print(f"Detected input format: {input_format}")
    
    if input_format == "unknown":
        raise ValueError(f"Could not determine format of input file: {input_path}")
    
    # Generate output path if not specified
    if output_path is None:
        suffix_map = {
            "nk": "_nk.txt",
            "dielectric": "_dielectric.txt", 
            "json": "_material.json"
        }
        suffix = suffix_map.get(output_format, f"_{output_format}.txt")
        output_path = input_path.parent / f"{input_path.stem}{suffix}"
    else:
        output_path = Path(output_path)
    
    print(f"Converting {input_path} -> {output_path}")
    print(f"Format: {input_format} -> {output_format}")
    
    # Load material
    try:
        if input_format == "nk":
            material = MaterialLoader.load_from_nk_file(
                input_path, 
                wavelength_unit=wavelength_unit
            )
        elif input_format == "dielectric":
            material = MaterialLoader.load_from_dielectric_file(
                input_path,
                energy_unit=energy_unit
            )
        elif input_format == "json":
            material = MaterialLoader.load_from_json(input_path)
        else:
            raise ValueError(f"Unsupported input format: {input_format}")
        
        print(f"Loaded material: {material}")
        
    except Exception as e:
        raise RuntimeError(f"Failed to load material: {e}")
    
    # Save in output format
    try:
        if output_format == "json":
            MaterialLoader.save_to_json(material, output_path)
        elif output_format == "dielectric":
            # Convert to dielectric format
            if input_format == "nk":
                MaterialLoader.convert_nk_to_dielectric(input_path, output_path, wavelength_unit)
            else:
                # Already have material object, save as dielectric
                save_material_as_dielectric(material, output_path)
        elif output_format == "nk":
            save_material_as_nk(material, output_path)
        else:
            raise ValueError(f"Unsupported output format: {output_format}")
        
        print(f"Conversion completed successfully!")
        print(f"Output saved to: {output_path}")
        
    except Exception as e:
        raise RuntimeError(f"Failed to save material: {e}")

def save_material_as_dielectric(material, output_path: Path) -> None:
    """Save material as dielectric function file."""
    import numpy as np
    
    # Get dielectric data
    optical_data = material.get_optical_constants(material.energy_eV)
    
    # Create output array
    data = np.column_stack([
        material.energy_eV,
        optical_data['eps_xx'].real, optical_data['eps_xx'].imag,
        optical_data['eps_yy'].real, optical_data['eps_yy'].imag,
        optical_data['eps_zz'].real, optical_data['eps_zz'].imag
    ])
    
    # Save with header
    header = "Energy(eV) eps_x_real eps_x_imag eps_y_real eps_y_imag eps_z_real eps_z_imag"
    np.savetxt(output_path, data, header=header, fmt="%.6e")

def save_material_as_nk(material, output_path: Path) -> None:
    """Save material as n,k file."""
    import numpy as np
    from ..utils.constants import Constants
    
    # Get optical data
    optical_data = material.get_optical_constants(material.energy_eV)
    
    # Convert energy to wavelength
    wavelength_nm = optical_data['wavelength_m'] * 1e9
    
    # Sort by wavelength (ascending)
    sort_idx = np.argsort(wavelength_nm)
    wavelength_sorted = wavelength_nm[sort_idx]
    n_sorted = optical_data['n_x'][sort_idx]
    k_sorted = optical_data['k_x'][sort_idx]
    
    # Create output array
    data = np.column_stack([wavelength_sorted, n_sorted, k_sorted])
    
    # Save with header
    header = "Wavelength(nm) n k"
    np.savetxt(output_path, data, header=header, fmt="%.6f")

def main():
    """Main command-line interface."""
    parser = argparse.ArgumentParser(
        description="Convert material files between different formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s material.txt                          # Auto-detect format, convert to dielectric
  %(prog)s material.txt -o material_eps.txt     # Specify output file
  %(prog)s material.txt -if nk -of json         # Convert n,k to JSON
  %(prog)s material.txt -of nk -wu nm           # Convert to n,k with nm wavelength
  %(prog)s material_eps.txt -if dielectric      # Explicitly specify input format
        """
    )
    
    parser.add_argument(
        "input", 
        help="Input material file"
    )
    
    parser.add_argument(
        "-o", "--output",
        help="Output file path (auto-generated if not specified)"
    )
    
    parser.add_argument(
        "-if", "--input-format",
        choices=["nk", "dielectric", "json"],
        help="Input file format (auto-detect if not specified)"
    )
    
    parser.add_argument(
        "-of", "--output-format", 
        choices=["nk", "dielectric", "json"],
        default="dielectric",
        help="Output file format (default: dielectric)"
    )
    
    parser.add_argument(
        "-wu", "--wavelength-unit",
        choices=["nm", "um", "m"],
        default="um",
        help="Wavelength unit for n,k files (default: um)"
    )
    
    parser.add_argument(
        "-eu", "--energy-unit",
        choices=["eV", "cm-1", "THz"],
        default="eV", 
        help="Energy unit for dielectric files (default: eV)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    try:
        convert_file(
            input_path=args.input,
            output_path=args.output,
            input_format=args.input_format,
            output_format=args.output_format,
            wavelength_unit=args.wavelength_unit,
            energy_unit=args.energy_unit
        )
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 