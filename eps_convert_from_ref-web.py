import numpy as np
import os
from tkinter import Tk, filedialog

def convert_dielectric(input_path):
    data = np.loadtxt(input_path)
    wavelength_um = data[:, 0]
    n = data[:, 1]
    k = data[:, 2]

    # Compute dielectric function
    eps_real = n**2 - k**2
    eps_imag = 2 * n * k

    # Convert wavelength (µm) to wavenumber (cm⁻¹)
    wavenumber = 1e4 / wavelength_um  # since 1 µm = 1e-4 cm

    # Stack output: wavenumber, eps1_real, eps1_im, eps2_real, ..., eps3_im
    output = np.column_stack([
        wavenumber,
        eps_real, eps_imag,
        eps_real, eps_imag,  # duplicate for y-axis
        eps_real, eps_imag   # duplicate for z-axis
    ])

    # Save to same folder
    folder, filename = os.path.split(input_path)
    base, _ = os.path.splitext(filename)
    output_path = os.path.join(folder, f"{base}_converted.txt")
    header = "wavenumber(cm^-1) epsx_real epsx_imag epsy_real epsy_imag epsz_real epsz_imag"
    np.savetxt(output_path, output, fmt="%.6f", header=header)

    print(f"Converted file saved to:\n{output_path}")

if __name__ == "__main__":
    Tk().withdraw()  # Hide the root window
    file_path = filedialog.askopenfilename(title="Select dielectric function TXT file",
                                           filetypes=[("Text files", "*.txt")])
    if file_path:
        convert_dielectric(file_path)
    else:
        print("No file selected.")