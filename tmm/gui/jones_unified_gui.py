"""
Unified Jones Calculator GUI

This module provides a comprehensive GUI interface for the unified Jones
matrix calculator, combining dispersion analysis, k-space cross-sections,
and 3D E-kx-ky analysis with adaptive resolution controls.

Author: Ding Xu (Physical Chemistry Researcher)
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading
import time

from ..core.layer import Layer, create_substrate, create_superstrate, create_film
from ..materials.builtin_materials import get_builtin_material
from ..calculations.jones_unified_calculator import (
    JonesUnifiedCalculator, CalculationParams, GridResolution
)
from ..calculations.jones_3d_analyzer import Jones3DAnalyzer


class JonesUnifiedGUI:
    """
    Unified GUI for Jones matrix calculations with adaptive resolution.
    """
    
    def __init__(self, root):
        self.root = root
        self.root.title("TMM - Unified Jones Calculator")
        self.root.geometry("1200x800")
        
        # Initialize variables
        self.layers = []
        self.results = None
        self.calculation_thread = None
        self.calculation_running = False
        
        # Create GUI
        self.create_widgets()
        self.setup_default_structure()
    
    def create_widgets(self):
        """Create the main GUI interface"""
        # Main notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Create tabs
        self.create_structure_tab()
        self.create_parameters_tab()
        self.create_calculation_tab()
        self.create_visualization_tab()
    
    def create_structure_tab(self):
        """Create layer structure configuration tab"""
        structure_frame = ttk.Frame(self.notebook)
        self.notebook.add(structure_frame, text="Layer Structure")
        
        # Layer table
        ttk.Label(structure_frame, text="Layer Structure", font=('Arial', 12, 'bold')).pack(pady=5)
        
        # Layer listbox with scrollbar
        list_frame = tk.Frame(structure_frame)
        list_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        scrollbar = tk.Scrollbar(list_frame)
        scrollbar.pack(side='right', fill='y')
        
        self.layer_listbox = tk.Listbox(list_frame, yscrollcommand=scrollbar.set, height=10)
        self.layer_listbox.pack(side='left', fill='both', expand=True)
        scrollbar.config(command=self.layer_listbox.yview)
        
        # Layer control buttons
        button_frame = tk.Frame(structure_frame)
        button_frame.pack(pady=5)
        
        tk.Button(button_frame, text="Add Layer", command=self.add_layer).pack(side='left', padx=2)
        tk.Button(button_frame, text="Edit Layer", command=self.edit_layer).pack(side='left', padx=2)
        tk.Button(button_frame, text="Delete Layer", command=self.delete_layer).pack(side='left', padx=2)
        tk.Button(button_frame, text="Load Structure", command=self.load_structure).pack(side='left', padx=2)
        tk.Button(button_frame, text="Save Structure", command=self.save_structure).pack(side='left', padx=2)
        
        # Predefined structures
        preset_frame = tk.Frame(structure_frame)
        preset_frame.pack(pady=5)
        
        ttk.Label(preset_frame, text="Preset Structures:").pack(side='left')
        tk.Button(preset_frame, text="Basic SiO2", command=self.load_basic_structure).pack(side='left', padx=2)
        tk.Button(preset_frame, text="Plasmonic", command=self.load_plasmonic_structure).pack(side='left', padx=2)
    
    def create_parameters_tab(self):
        """Create calculation parameters tab"""
        params_frame = ttk.Frame(self.notebook)
        self.notebook.add(params_frame, text="Parameters")
        
        # Create parameter sections
        self.create_energy_params(params_frame)
        self.create_kspace_params(params_frame)
        self.create_polarization_params(params_frame)
        self.create_resolution_params(params_frame)
    
    def create_energy_params(self, parent):
        """Create energy parameter controls"""
        energy_frame = ttk.LabelFrame(parent, text="Energy Parameters", padding=10)
        energy_frame.pack(fill='x', padx=5, pady=5)
        
        # Energy range
        row1 = tk.Frame(energy_frame)
        row1.pack(fill='x', pady=2)
        
        tk.Label(row1, text="Energy Start (eV):").pack(side='left')
        self.energy_start = tk.DoubleVar(value=1.0)
        tk.Entry(row1, textvariable=self.energy_start, width=10).pack(side='left', padx=5)
        
        tk.Label(row1, text="Energy Stop (eV):").pack(side='left', padx=(20,0))
        self.energy_stop = tk.DoubleVar(value=3.0)
        tk.Entry(row1, textvariable=self.energy_stop, width=10).pack(side='left', padx=5)
        
        tk.Label(row1, text="Points:").pack(side='left', padx=(20,0))
        self.energy_points = tk.IntVar(value=101)
        tk.Entry(row1, textvariable=self.energy_points, width=8).pack(side='left', padx=5)
    
    def create_kspace_params(self, parent):
        """Create k-space parameter controls"""
        kspace_frame = ttk.LabelFrame(parent, text="k-Space Parameters", padding=10)
        kspace_frame.pack(fill='x', padx=5, pady=5)
        
        # kx parameters
        row1 = tk.Frame(kspace_frame)
        row1.pack(fill='x', pady=2)
        
        tk.Label(row1, text="kx Start (μm⁻¹):").pack(side='left')
        self.kx_start = tk.DoubleVar(value=0.0)
        tk.Entry(row1, textvariable=self.kx_start, width=10).pack(side='left', padx=5)
        
        tk.Label(row1, text="kx Stop (μm⁻¹):").pack(side='left', padx=(20,0))
        self.kx_stop = tk.DoubleVar(value=25.0)
        tk.Entry(row1, textvariable=self.kx_stop, width=10).pack(side='left', padx=5)
        
        tk.Label(row1, text="Points:").pack(side='left', padx=(20,0))
        self.kx_points = tk.IntVar(value=101)
        tk.Entry(row1, textvariable=self.kx_points, width=8).pack(side='left', padx=5)
        
        # ky parameters
        row2 = tk.Frame(kspace_frame)
        row2.pack(fill='x', pady=2)
        
        tk.Label(row2, text="ky Start (μm⁻¹):").pack(side='left')
        self.ky_start = tk.DoubleVar(value=0.0)
        tk.Entry(row2, textvariable=self.ky_start, width=10).pack(side='left', padx=5)
        
        tk.Label(row2, text="ky Stop (μm⁻¹):").pack(side='left', padx=(20,0))
        self.ky_stop = tk.DoubleVar(value=25.0)
        tk.Entry(row2, textvariable=self.ky_stop, width=10).pack(side='left', padx=5)
        
        tk.Label(row2, text="Points:").pack(side='left', padx=(20,0))
        self.ky_points = tk.IntVar(value=101)
        tk.Entry(row2, textvariable=self.ky_points, width=8).pack(side='left', padx=5)
    
    def create_polarization_params(self, parent):
        """Create polarization parameter controls"""
        pol_frame = ttk.LabelFrame(parent, text="Polarization Parameters", padding=10)
        pol_frame.pack(fill='x', padx=5, pady=5)
        
        row = tk.Frame(pol_frame)
        row.pack(fill='x', pady=2)
        
        tk.Label(row, text="|Es| (s-amplitude):").pack(side='left')
        self.es_amplitude = tk.DoubleVar(value=1.0)
        tk.Entry(row, textvariable=self.es_amplitude, width=8).pack(side='left', padx=5)
        
        tk.Label(row, text="|Ep| (p-amplitude):").pack(side='left', padx=(20,0))
        self.ep_amplitude = tk.DoubleVar(value=0.0)
        tk.Entry(row, textvariable=self.ep_amplitude, width=8).pack(side='left', padx=5)
        
        tk.Label(row, text="Phase Diff (°):").pack(side='left', padx=(20,0))
        self.phase_diff = tk.DoubleVar(value=0.0)
        tk.Entry(row, textvariable=self.phase_diff, width=8).pack(side='left', padx=5)
    
    def create_resolution_params(self, parent):
        """Create resolution parameter controls"""
        res_frame = ttk.LabelFrame(parent, text="Resolution & Performance", padding=10)
        res_frame.pack(fill='x', padx=5, pady=5)
        
        row1 = tk.Frame(res_frame)
        row1.pack(fill='x', pady=2)
        
        tk.Label(row1, text="Grid Resolution:").pack(side='left')
        self.resolution = tk.StringVar(value="medium")
        resolution_combo = ttk.Combobox(row1, textvariable=self.resolution, 
                                      values=["coarse", "medium", "fine", "custom"], 
                                      state="readonly", width=10)
        resolution_combo.pack(side='left', padx=5)
        
        tk.Label(row1, text="Max Memory (GB):").pack(side='left', padx=(20,0))
        self.max_memory = tk.DoubleVar(value=8.0)
        tk.Entry(row1, textvariable=self.max_memory, width=8).pack(side='left', padx=5)
        
        row2 = tk.Frame(res_frame)
        row2.pack(fill='x', pady=2)
        
        self.use_chunking = tk.BooleanVar(value=True)
        tk.Checkbutton(row2, text="Use chunking for large calculations", 
                      variable=self.use_chunking).pack(side='left')
        
        # Memory estimation
        tk.Button(row2, text="Estimate Memory", 
                 command=self.estimate_memory).pack(side='right', padx=5)
    
    def create_calculation_tab(self):
        """Create calculation control tab"""
        calc_frame = ttk.Frame(self.notebook)
        self.notebook.add(calc_frame, text="Calculation")
        
        # Calculation type selection
        type_frame = ttk.LabelFrame(calc_frame, text="Calculation Type", padding=10)
        type_frame.pack(fill='x', padx=5, pady=5)
        
        self.calc_type = tk.StringVar(value="3d_full")
        tk.Radiobutton(type_frame, text="3D E-kx-ky Full Analysis", 
                      variable=self.calc_type, value="3d_full").pack(anchor='w')
        tk.Radiobutton(type_frame, text="2D kx-ky Cross-section", 
                      variable=self.calc_type, value="2d_kspace").pack(anchor='w')
        tk.Radiobutton(type_frame, text="1D Dispersion Line", 
                      variable=self.calc_type, value="1d_dispersion").pack(anchor='w')
        
        # Specific parameters for each type
        spec_frame = ttk.LabelFrame(calc_frame, text="Specific Parameters", padding=10)
        spec_frame.pack(fill='x', padx=5, pady=5)
        
        # 2D slice energy
        row1 = tk.Frame(spec_frame)
        row1.pack(fill='x', pady=2)
        tk.Label(row1, text="2D Slice Energy (eV):").pack(side='left')
        self.slice_energy = tk.DoubleVar(value=2.0)
        tk.Entry(row1, textvariable=self.slice_energy, width=10).pack(side='left', padx=5)
        
        # 1D dispersion parameters
        row2 = tk.Frame(spec_frame)
        row2.pack(fill='x', pady=2)
        tk.Label(row2, text="1D Fixed Axis:").pack(side='left')
        self.fixed_axis = tk.StringVar(value="kx")
        ttk.Combobox(row2, textvariable=self.fixed_axis, values=["kx", "ky"], 
                    state="readonly", width=8).pack(side='left', padx=5)
        
        tk.Label(row2, text="Fixed Value (μm⁻¹):").pack(side='left', padx=(20,0))
        self.fixed_value = tk.DoubleVar(value=0.0)
        tk.Entry(row2, textvariable=self.fixed_value, width=10).pack(side='left', padx=5)
        
        # Control buttons
        control_frame = tk.Frame(calc_frame)
        control_frame.pack(pady=20)
        
        self.calc_button = tk.Button(control_frame, text="Start Calculation", 
                                   command=self.start_calculation, font=('Arial', 12, 'bold'))
        self.calc_button.pack(side='left', padx=5)
        
        self.stop_button = tk.Button(control_frame, text="Stop Calculation", 
                                   command=self.stop_calculation, state='disabled')
        self.stop_button.pack(side='left', padx=5)
        
        # Progress bar
        self.progress = ttk.Progressbar(calc_frame, mode='indeterminate')
        self.progress.pack(fill='x', padx=5, pady=5)
        
        # Status text
        self.status_text = tk.Text(calc_frame, height=8, wrap='word')
        self.status_text.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Scrollbar for status
        status_scroll = tk.Scrollbar(self.status_text)
        status_scroll.pack(side='right', fill='y')
        self.status_text.config(yscrollcommand=status_scroll.set)
        self.status_text.config(command=self.status_text.yview)
    
    def create_visualization_tab(self):
        """Create visualization tab"""
        vis_frame = ttk.Frame(self.notebook)
        self.notebook.add(vis_frame, text="Visualization")
        
        # Placeholder for visualization controls
        tk.Label(vis_frame, text="Visualization controls will be added here", 
                font=('Arial', 12)).pack(pady=20)
        
        # Quick plot buttons
        plot_frame = tk.Frame(vis_frame)
        plot_frame.pack(pady=10)
        
        tk.Button(plot_frame, text="Plot Results", command=self.plot_results).pack(side='left', padx=5)
        tk.Button(plot_frame, text="Save Data", command=self.save_data).pack(side='left', padx=5)
        tk.Button(plot_frame, text="Export Plots", command=self.export_plots).pack(side='left', padx=5)
    
    def setup_default_structure(self):
        """Set up default layer structure"""
        self.load_basic_structure()
    
    def load_basic_structure(self):
        """Load basic SiO2 structure"""
        try:
            air = get_builtin_material("Air")
            sio2 = get_builtin_material("SiO2")
            
            self.layers = [
                create_superstrate(air, "Air (top)"),
                create_film(sio2, 500, "SiO2 500nm"),
                create_substrate(air, "Air (bottom)")
            ]
            
            self.update_layer_display()
            self.log_message("Loaded basic SiO2 structure")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load basic structure: {str(e)}")
    
    def load_plasmonic_structure(self):
        """Load plasmonic structure"""
        try:
            air = get_builtin_material("Air")
            gold = get_builtin_material("Gold")
            sio2 = get_builtin_material("SiO2")
            
            self.layers = [
                create_superstrate(air, "Air (top)"),
                create_film(gold, 50, "Gold 50nm"),
                create_substrate(sio2, "SiO2 (substrate)")
            ]
            
            self.update_layer_display()
            self.log_message("Loaded plasmonic structure")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load plasmonic structure: {str(e)}")
    
    def update_layer_display(self):
        """Update the layer display"""
        self.layer_listbox.delete(0, tk.END)
        for i, layer in enumerate(self.layers):
            if layer.is_semi_infinite:
                thickness_str = "∞"
            else:
                thickness_str = f"{layer.thickness_nm:.1f} nm"
            
            display_str = f"{i+1}. {layer.name} ({thickness_str})"
            self.layer_listbox.insert(tk.END, display_str)
    
    def add_layer(self):
        """Add a new layer"""
        # This would open a dialog to add a new layer
        messagebox.showinfo("Add Layer", "Layer addition dialog will be implemented")
    
    def edit_layer(self):
        """Edit selected layer"""
        messagebox.showinfo("Edit Layer", "Layer editing dialog will be implemented")
    
    def delete_layer(self):
        """Delete selected layer"""
        messagebox.showinfo("Delete Layer", "Layer deletion will be implemented")
    
    def load_structure(self):
        """Load structure from file"""
        messagebox.showinfo("Load Structure", "Structure loading will be implemented")
    
    def save_structure(self):
        """Save structure to file"""
        messagebox.showinfo("Save Structure", "Structure saving will be implemented")
    
    def estimate_memory(self):
        """Estimate memory usage for current parameters"""
        try:
            params = self.get_calculation_params()
            
            # Create temporary calculator for estimation
            calculator = JonesUnifiedCalculator(self.layers)
            memory_mb = calculator.estimate_memory_usage(params)
            
            messagebox.showinfo("Memory Estimate", 
                              f"Estimated memory usage: {memory_mb:.1f} MB\n"
                              f"Grid size: {params.energy_points} × {params.kx_points} × {params.ky_points}\n"
                              f"Total points: {params.energy_points * params.kx_points * params.ky_points}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to estimate memory: {str(e)}")
    
    def get_calculation_params(self):
        """Get calculation parameters from GUI"""
        return CalculationParams(
            energy_start_eV=self.energy_start.get(),
            energy_stop_eV=self.energy_stop.get(),
            energy_points=self.energy_points.get(),
            kx_start_um=self.kx_start.get(),
            kx_stop_um=self.kx_stop.get(),
            kx_points=self.kx_points.get(),
            ky_start_um=self.ky_start.get(),
            ky_stop_um=self.ky_stop.get(),
            ky_points=self.ky_points.get(),
            es_amplitude=self.es_amplitude.get(),
            ep_amplitude=self.ep_amplitude.get(),
            phase_diff_deg=self.phase_diff.get(),
            resolution=GridResolution(self.resolution.get()),
            max_memory_gb=self.max_memory.get(),
            use_chunking=self.use_chunking.get()
        )
    
    def start_calculation(self):
        """Start calculation in background thread"""
        if self.calculation_running:
            messagebox.showwarning("Warning", "Calculation already running")
            return
        
        if not self.layers:
            messagebox.showerror("Error", "No layer structure defined")
            return
        
        try:
            params = self.get_calculation_params()
            calc_type = self.calc_type.get()
            
            # Update UI
            self.calc_button.config(state='disabled')
            self.stop_button.config(state='normal')
            self.progress.start()
            self.calculation_running = True
            
            # Start calculation thread
            self.calculation_thread = threading.Thread(
                target=self.run_calculation, 
                args=(params, calc_type)
            )
            self.calculation_thread.start()
            
        except Exception as e:
            self.log_message(f"Error starting calculation: {str(e)}")
            self.reset_calculation_ui()
    
    def run_calculation(self, params, calc_type):
        """Run calculation in background thread"""
        try:
            self.log_message(f"Starting {calc_type} calculation...")
            
            # Create analyzer
            analyzer = Jones3DAnalyzer(self.layers)
            
            if calc_type == "3d_full":
                self.results = analyzer.calculate_3d_full(params)
            elif calc_type == "2d_kspace":
                # First calculate 3D with single energy
                params_2d = params
                params_2d.energy_start_eV = self.slice_energy.get()
                params_2d.energy_stop_eV = self.slice_energy.get()
                params_2d.energy_points = 1
                self.results = analyzer.calculate_3d_full(params_2d)
            elif calc_type == "1d_dispersion":
                # This would be implemented for 1D dispersion
                self.log_message("1D dispersion calculation not yet implemented")
                return
            
            self.log_message(f"Calculation completed successfully!")
            self.log_message(f"Computation time: {self.results.computation_time:.2f} seconds")
            
            # Update UI in main thread
            self.root.after(0, self.calculation_finished)
            
        except Exception as e:
            self.log_message(f"Calculation failed: {str(e)}")
            self.root.after(0, self.reset_calculation_ui)
    
    def calculation_finished(self):
        """Handle calculation completion"""
        self.reset_calculation_ui()
        self.notebook.select(3)  # Switch to visualization tab
        messagebox.showinfo("Calculation Complete", 
                          f"Calculation completed successfully!\n"
                          f"Time: {self.results.computation_time:.2f} seconds")
    
    def stop_calculation(self):
        """Stop running calculation"""
        self.calculation_running = False
        self.reset_calculation_ui()
        self.log_message("Calculation stopped by user")
    
    def reset_calculation_ui(self):
        """Reset calculation UI elements"""
        self.calc_button.config(state='normal')
        self.stop_button.config(state='disabled')
        self.progress.stop()
        self.calculation_running = False
    
    def log_message(self, message):
        """Log message to status text"""
        timestamp = time.strftime("%H:%M:%S")
        self.status_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.status_text.see(tk.END)
        self.root.update()
    
    def plot_results(self):
        """Plot calculation results"""
        if self.results is None:
            messagebox.showwarning("Warning", "No results to plot")
            return
        
        try:
            # Create a simple plot
            fig, axes = plt.subplots(2, 4, figsize=(16, 8))
            
            # Plot all 8 quantities
            quantities = ['Rs', 'Rp', 'S0', 'S1', 'S2', 'S3', 'phi', 'chi']
            data_arrays = [self.results.Rs, self.results.Rp, self.results.S0, 
                          self.results.S1, self.results.S2, self.results.S3,
                          self.results.phi, self.results.chi]
            
            for i, (ax, title, data) in enumerate(zip(axes.flat, quantities, data_arrays)):
                if data.ndim == 3:
                    # 3D data - plot middle energy slice
                    mid_energy = data.shape[0] // 2
                    im = ax.pcolormesh(self.results.kx_um, self.results.ky_um, 
                                     data[mid_energy, :, :], shading='auto')
                    ax.set_title(f'{title} @ {self.results.energy_eV[mid_energy]:.2f} eV')
                    ax.set_xlabel('kx (μm⁻¹)')
                    ax.set_ylabel('ky (μm⁻¹)')
                elif data.ndim == 2:
                    # 2D data
                    im = ax.pcolormesh(self.results.kx_um, self.results.ky_um, 
                                     data, shading='auto')
                    ax.set_title(title)
                    ax.set_xlabel('kx (μm⁻¹)')
                    ax.set_ylabel('ky (μm⁻¹)')
                
                plt.colorbar(im, ax=ax)
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to plot results: {str(e)}")
    
    def save_data(self):
        """Save calculation data"""
        if self.results is None:
            messagebox.showwarning("Warning", "No results to save")
            return
        
        messagebox.showinfo("Save Data", "Data saving will be implemented")
    
    def export_plots(self):
        """Export plots to files"""
        if self.results is None:
            messagebox.showwarning("Warning", "No results to export")
            return
        
        messagebox.showinfo("Export Plots", "Plot export will be implemented")


def main():
    """Main function to run the GUI"""
    root = tk.Tk()
    app = JonesUnifiedGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main() 