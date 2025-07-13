import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, filedialog
import os
import time
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Physical constants
c = 3e8                    # speed of light (m/s)
mu0 = 4e-7 * np.pi         # vacuum permeability
eps0 = 1/(mu0 * c**2)      # vacuum permittivity

# Material database (name -> filepath)
material_files = {
    'Air': None,
    'SiO2': None
}  # type: dict[str, str | None]

# --- Material routines ---
def get_material_eps(material, w_cm):
    """Get material dielectric tensor for given frequencies"""
    if material == 'Air':
        eps = np.ones((len(w_cm), 3), dtype=complex)
    elif material == 'SiO2':
        n = 1.5
        eps = np.full((len(w_cm), 3), n**2, dtype=complex)
    elif material in material_files and material_files[material] is not None:
        filepath = material_files[material]
        assert filepath is not None  # Type guard for linter
        data = np.loadtxt(filepath)
        w_data = data[:, 0]
        eps = np.zeros((len(w_cm), 3), dtype=complex)
        for i, comp in enumerate([1, 3, 5]):
            re = np.interp(w_cm, w_data, data[:, comp])
            im = np.interp(w_cm, w_data, data[:, comp + 1])
            eps[:, i] = re + 1j * im
    else:
        raise ValueError(f"Unknown material or missing file: {material}")
    return eps.reshape(-1, 3, 1)

def build_layers(layer_info, w_cm):
    """Build layer structure with dielectric tensors"""
    layers = []
    for name, thickness in layer_info:
        layers.append({
            'd': thickness * 1e-9,
            'eps': get_material_eps(name, w_cm)
        })
    return layers

# --- Jones calculus routines ---
def local_basis_rotation(kx, ky):
    """Rotation matrix from lab frame to local s-p frame"""
    k_par = np.sqrt(kx**2 + ky**2)
    if k_par == 0:
        return np.eye(3)
    s_hat = np.array([-ky, kx, 0]) / k_par
    z_hat = np.array([0, 0, 1])
    p_hat = np.cross(z_hat, s_hat)
    return np.column_stack([s_hat, p_hat, z_hat])

def rotate_eps_full(eps_vec, kx, ky):
    """Rotate dielectric tensor to local s-p coordinates"""
    R = local_basis_rotation(kx, ky)
    eps_mat = np.diag(eps_vec)
    return R.T @ eps_mat @ R

def interface_jones(eps1, eps2, k_par, w):
    """Calculate Jones reflection matrix at interface"""
    # Normal components for s and p polarizations
    kz1s = np.sqrt(eps1[0,0]*(w/c)**2 - k_par**2 + 0j)
    kz2s = np.sqrt(eps2[0,0]*(w/c)**2 - k_par**2 + 0j)
    kz1p = np.sqrt(eps1[1,1]*(w/c)**2 - (eps1[1,1]/eps1[2,2])*k_par**2 + 0j)
    kz2p = np.sqrt(eps2[1,1]*(w/c)**2 - (eps2[1,1]/eps2[2,2])*k_par**2 + 0j)
    
    # Admittances
    Ys1 = kz1s/(mu0*w);  Yp1 = eps0*eps1[2,2]*w/kz1p
    Ys2 = kz2s/(mu0*w);  Yp2 = eps0*eps2[2,2]*w/kz2p
    Y1 = np.diag([Ys1, Yp1]); Y2 = np.diag([Ys2, Yp2])
    
    # Jones reflection matrix R = (Y2-Y1) @ inv(Y2+Y1)
    return np.linalg.solve(Y2+Y1, Y2-Y1)

def gamma_jones(rj, G_prev, P):
    """Cascade Jones matrices through layer"""
    # P = diag([e^{2iφ_s}, e^{2iφ_p}])
    num = rj + P @ G_prev @ P
    den = np.eye(2) + rj @ P @ G_prev @ P
    return np.linalg.solve(den, num)

# --- Main 3D calculation ---
def calc_3d_reflectance(layer_info, E_start, E_stop, dE, 
                       kx_start, kx_stop, dkx, ky_start, ky_stop, dky,
                       a, b, delta_phi_deg, progress_callback=None):
    """
    Calculate 3D reflectance and Stokes parameters in E-kx-ky space
    
    Parameters:
    - layer_info: List of (material, thickness) tuples
    - E_start, E_stop, dE: Energy range and step (eV)
    - kx_start, kx_stop, dkx: kx range and step (μm⁻¹)
    - ky_start, ky_stop, dky: ky range and step (μm⁻¹)
    - a, b: Input field amplitudes |Es|, |Ep|
    - delta_phi_deg: Phase difference in degrees
    - progress_callback: Function to call with progress updates
    """
    
    # Create grids
    E_vals = np.arange(E_start, E_stop + dE, dE)
    kx_vals = np.arange(kx_start, kx_stop + dkx, dkx)
    ky_vals = np.arange(ky_start, ky_stop + dky, dky)
    
    nE, nKx, nKy = len(E_vals), len(kx_vals), len(ky_vals)
    
    # Initialize output arrays
    Rs = np.zeros((nE, nKx, nKy))
    Rp = np.zeros((nE, nKx, nKy))
    S0 = np.zeros((nE, nKx, nKy))
    S1 = np.zeros((nE, nKx, nKy))
    S2 = np.zeros((nE, nKx, nKy))
    S3 = np.zeros((nE, nKx, nKy))
    phi = np.zeros((nE, nKx, nKy))
    chi = np.zeros((nE, nKx, nKy))
    
    # Input field in lab frame
    delta = np.deg2rad(delta_phi_deg)
    Ein_lab = np.array([-b*np.exp(1j*delta), a, 0])
    
    total_iterations = nE * nKx * nKy
    current_iteration = 0
    
    # Main calculation loop
    for i, E in enumerate(E_vals):
        # Convert energy to frequency
        w_spec = E * 8065.54429  # eV to cm⁻¹
        w = 2*np.pi * c * w_spec * 100  # angular frequency (rad/s)
        
        # Build layers for this energy
        layers = build_layers(layer_info, np.array([w_spec]))
        
        for j, kx in enumerate(kx_vals):
            for k, ky in enumerate(ky_vals):
                # Convert to m⁻¹
                kx_m = kx * 1e6
                ky_m = ky * 1e6
                k_par = np.sqrt(kx_m**2 + ky_m**2)
                
                # Calculate Jones reflection matrix by cascading interfaces
                G = None
                for iface in reversed(range(len(layers)-1)):
                    # Get rotated dielectric tensors
                    eps1 = rotate_eps_full(layers[iface]['eps'][0,:,0], kx_m, ky_m)
                    eps2 = rotate_eps_full(layers[iface+1]['eps'][0,:,0], kx_m, ky_m)
                    
                    # Interface reflection matrix
                    rj = interface_jones(eps1, eps2, k_par, w)
                    
                    if iface == len(layers) - 2:  # Last interface
                        G = rj
                    else:  # Cascade with previous layers
                        # Propagation phases
                        kz_s = np.sqrt(eps2[0,0]*(w/c)**2 - k_par**2 + 0j)
                        kz_p = np.sqrt(eps2[1,1]*(w/c)**2 - (eps2[1,1]/eps2[2,2])*k_par**2 + 0j)
                        phi_s = kz_s * layers[iface+1]['d']
                        phi_p = kz_p * layers[iface+1]['d']
                        phi_prop = np.diag([np.exp(1j*phi_s), np.exp(1j*phi_p)])
                        
                        # Cascade
                        G = gamma_jones(rj, G, phi_prop)
                
                # Handle case where G might be None (shouldn't happen with proper layer structure)
                if G is None:
                    G = np.zeros((2, 2), dtype=complex)
                
                # Apply Jones matrix to input field
                R3 = local_basis_rotation(kx_m, ky_m)
                Ein_loc = R3.T @ Ein_lab
                Eout_loc = G @ Ein_loc[:2]
                Eref_lab = R3 @ np.array([Eout_loc[0], Eout_loc[1], 0])
                
                # Extract field components
                Ex, Ey = Eref_lab[0], Eref_lab[1]
                
                # Calculate Stokes parameters
                S0_val = abs(Ex)**2 + abs(Ey)**2
                if S0_val > 0:
                    S1_val = (abs(Ex)**2 - abs(Ey)**2) / S0_val
                    S2_val = 2*np.real(Ex*np.conj(Ey)) / S0_val
                    S3_val = 2*np.imag(Ex*np.conj(Ey)) / S0_val
                    phi_val = 0.5*np.arctan2(S2_val, S1_val)
                    chi_val = 0.5*np.arcsin(np.clip(S3_val, -1, 1))
                else:
                    S1_val = S2_val = S3_val = phi_val = chi_val = 0
                
                # Store results
                Rs[i,j,k] = abs(G[0,0])**2
                Rp[i,j,k] = abs(G[1,1])**2
                S0[i,j,k] = S0_val
                S1[i,j,k] = S1_val
                S2[i,j,k] = S2_val
                S3[i,j,k] = S3_val
                phi[i,j,k] = phi_val
                chi[i,j,k] = chi_val
                
                # Update progress
                current_iteration += 1
                if progress_callback and current_iteration % 100 == 0:
                    progress = current_iteration / total_iterations
                    progress_callback(progress)
    
    return E_vals, kx_vals, ky_vals, Rs, Rp, S0, S1, S2, S3, phi, chi

# --- GUI Implementation ---
class TMM3DGui:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title('3D TMM Calculator - E-kx-ky')
        self.root.geometry('800x700')
        
        # Initialize layer structure
        self.layers = [('Air', 1e6), ('SiO2', 500), ('Air', 1e6)]
        self.materials = list(material_files.keys())
        
        # Results storage
        self.results = None
        
        self.setup_gui()
        
    def setup_gui(self):
        # Main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Parameter input frame
        param_frame = ttk.LabelFrame(main_frame, text="Calculation Parameters")
        param_frame.pack(fill='x', pady=(0, 10))
        
        # Energy parameters
        energy_frame = ttk.Frame(param_frame)
        energy_frame.pack(fill='x', pady=5)
        
        ttk.Label(energy_frame, text='Energy (eV):').grid(row=0, column=0, padx=5)
        ttk.Label(energy_frame, text='Start:').grid(row=0, column=1, padx=5)
        self.E_start = tk.DoubleVar(value=1.0)
        ttk.Entry(energy_frame, textvariable=self.E_start, width=8).grid(row=0, column=2, padx=5)
        
        ttk.Label(energy_frame, text='Stop:').grid(row=0, column=3, padx=5)
        self.E_stop = tk.DoubleVar(value=3.0)
        ttk.Entry(energy_frame, textvariable=self.E_stop, width=8).grid(row=0, column=4, padx=5)
        
        ttk.Label(energy_frame, text='Step:').grid(row=0, column=5, padx=5)
        self.dE = tk.DoubleVar(value=0.1)
        ttk.Entry(energy_frame, textvariable=self.dE, width=8).grid(row=0, column=6, padx=5)
        
        # kx parameters
        kx_frame = ttk.Frame(param_frame)
        kx_frame.pack(fill='x', pady=5)
        
        ttk.Label(kx_frame, text='kx (μm⁻¹):').grid(row=0, column=0, padx=5)
        ttk.Label(kx_frame, text='Start:').grid(row=0, column=1, padx=5)
        self.kx_start = tk.DoubleVar(value=-10.0)
        ttk.Entry(kx_frame, textvariable=self.kx_start, width=8).grid(row=0, column=2, padx=5)
        
        ttk.Label(kx_frame, text='Stop:').grid(row=0, column=3, padx=5)
        self.kx_stop = tk.DoubleVar(value=10.0)
        ttk.Entry(kx_frame, textvariable=self.kx_stop, width=8).grid(row=0, column=4, padx=5)
        
        ttk.Label(kx_frame, text='Step:').grid(row=0, column=5, padx=5)
        self.dkx = tk.DoubleVar(value=1.0)
        ttk.Entry(kx_frame, textvariable=self.dkx, width=8).grid(row=0, column=6, padx=5)
        
        # ky parameters
        ky_frame = ttk.Frame(param_frame)
        ky_frame.pack(fill='x', pady=5)
        
        ttk.Label(ky_frame, text='ky (μm⁻¹):').grid(row=0, column=0, padx=5)
        ttk.Label(ky_frame, text='Start:').grid(row=0, column=1, padx=5)
        self.ky_start = tk.DoubleVar(value=-10.0)
        ttk.Entry(ky_frame, textvariable=self.ky_start, width=8).grid(row=0, column=2, padx=5)
        
        ttk.Label(ky_frame, text='Stop:').grid(row=0, column=3, padx=5)
        self.ky_stop = tk.DoubleVar(value=10.0)
        ttk.Entry(ky_frame, textvariable=self.ky_stop, width=8).grid(row=0, column=4, padx=5)
        
        ttk.Label(ky_frame, text='Step:').grid(row=0, column=5, padx=5)
        self.dky = tk.DoubleVar(value=1.0)
        ttk.Entry(ky_frame, textvariable=self.dky, width=8).grid(row=0, column=6, padx=5)
        
        # Polarization parameters
        pol_frame = ttk.Frame(param_frame)
        pol_frame.pack(fill='x', pady=5)
        
        ttk.Label(pol_frame, text='Input Polarization:').grid(row=0, column=0, padx=5)
        ttk.Label(pol_frame, text='|Es|:').grid(row=0, column=1, padx=5)
        self.Es_amp = tk.DoubleVar(value=1.0)
        ttk.Entry(pol_frame, textvariable=self.Es_amp, width=8).grid(row=0, column=2, padx=5)
        
        ttk.Label(pol_frame, text='|Ep|:').grid(row=0, column=3, padx=5)
        self.Ep_amp = tk.DoubleVar(value=0.0)
        ttk.Entry(pol_frame, textvariable=self.Ep_amp, width=8).grid(row=0, column=4, padx=5)
        
        ttk.Label(pol_frame, text='Δφ (°):').grid(row=0, column=5, padx=5)
        self.delta_phi = tk.DoubleVar(value=0.0)
        ttk.Entry(pol_frame, textvariable=self.delta_phi, width=8).grid(row=0, column=6, padx=5)
        
        # Grid size preset buttons
        preset_frame = ttk.Frame(param_frame)
        preset_frame.pack(fill='x', pady=5)
        
        ttk.Label(preset_frame, text='Grid Presets:').pack(side='left', padx=5)
        ttk.Button(preset_frame, text='Coarse (Test)', command=self.set_coarse_grid).pack(side='left', padx=5)
        ttk.Button(preset_frame, text='Medium', command=self.set_medium_grid).pack(side='left', padx=5)
        ttk.Button(preset_frame, text='Fine (Final)', command=self.set_fine_grid).pack(side='left', padx=5)
        
        # Layer management frame
        layer_frame = ttk.LabelFrame(main_frame, text="Layer Structure")
        layer_frame.pack(fill='both', expand=True, pady=(0, 10))
        
        # Layer table
        self.layer_table = ttk.Treeview(layer_frame, columns=('Material', 'Thickness'), show='headings', height=6)
        self.layer_table.heading('Material', text='Material')
        self.layer_table.heading('Thickness', text='Thickness (nm)')
        self.layer_table.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Initialize table
        self.update_layer_table()
        
        # Layer control buttons
        layer_btn_frame = ttk.Frame(layer_frame)
        layer_btn_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(layer_btn_frame, text='Add Layer', command=self.add_layer).pack(side='left', padx=5)
        ttk.Button(layer_btn_frame, text='Edit Layer', command=self.edit_layer).pack(side='left', padx=5)
        ttk.Button(layer_btn_frame, text='Delete Layer', command=self.delete_layer).pack(side='left', padx=5)
        ttk.Button(layer_btn_frame, text='Load Material File', command=self.load_material_file).pack(side='left', padx=5)
        
        # Control buttons frame
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(fill='x', pady=5)
        
        ttk.Button(control_frame, text='Calculate 3D', command=self.calculate_3d).pack(side='left', padx=5)
        ttk.Button(control_frame, text='Visualize Results', command=self.visualize_results).pack(side='left', padx=5)
        ttk.Button(control_frame, text='Save Results', command=self.save_results).pack(side='left', padx=5)
        
        # Progress bar
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(main_frame, variable=self.progress_var, maximum=100)
        self.progress_bar.pack(fill='x', pady=5)
        
        # Status label
        self.status_label = ttk.Label(main_frame, text="Ready")
        self.status_label.pack(fill='x', pady=5)
        
    def set_coarse_grid(self):
        """Set coarse grid for testing"""
        self.dE.set(0.5)
        self.dkx.set(2.0)
        self.dky.set(2.0)
        
    def set_medium_grid(self):
        """Set medium grid"""
        self.dE.set(0.2)
        self.dkx.set(1.0)
        self.dky.set(1.0)
        
    def set_fine_grid(self):
        """Set fine grid for final results"""
        self.dE.set(0.05)
        self.dkx.set(0.5)
        self.dky.set(0.5)
        
    def update_layer_table(self):
        """Update the layer table display"""
        for item in self.layer_table.get_children():
            self.layer_table.delete(item)
        for material, thickness in self.layers:
            self.layer_table.insert('', 'end', values=(material, thickness))
            
    def add_layer(self):
        """Add a new layer"""
        dialog = tk.Toplevel(self.root)
        dialog.title('Add Layer')
        dialog.geometry('300x150')
        
        ttk.Label(dialog, text='Material:').grid(row=0, column=0, padx=5, pady=5)
        material_var = tk.StringVar(value=self.materials[0])
        material_combo = ttk.Combobox(dialog, textvariable=material_var, values=self.materials, state='readonly')
        material_combo.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(dialog, text='Thickness (nm):').grid(row=1, column=0, padx=5, pady=5)
        thickness_var = tk.DoubleVar(value=100.0)
        thickness_entry = ttk.Entry(dialog, textvariable=thickness_var)
        thickness_entry.grid(row=1, column=1, padx=5, pady=5)
        
        def add_layer_ok():
            self.layers.append((material_var.get(), thickness_var.get()))
            self.update_layer_table()
            dialog.destroy()
            
        ttk.Button(dialog, text='Add', command=add_layer_ok).grid(row=2, column=0, columnspan=2, pady=10)
        
    def edit_layer(self):
        """Edit selected layer"""
        selected = self.layer_table.selection()
        if not selected:
            messagebox.showwarning("Warning", "Please select a layer to edit")
            return
            
        index = self.layer_table.index(selected[0])
        old_material, old_thickness = self.layers[index]
        
        dialog = tk.Toplevel(self.root)
        dialog.title('Edit Layer')
        dialog.geometry('300x150')
        
        ttk.Label(dialog, text='Material:').grid(row=0, column=0, padx=5, pady=5)
        material_var = tk.StringVar(value=old_material)
        material_combo = ttk.Combobox(dialog, textvariable=material_var, values=self.materials, state='readonly')
        material_combo.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(dialog, text='Thickness (nm):').grid(row=1, column=0, padx=5, pady=5)
        thickness_var = tk.DoubleVar(value=old_thickness)
        thickness_entry = ttk.Entry(dialog, textvariable=thickness_var)
        thickness_entry.grid(row=1, column=1, padx=5, pady=5)
        
        def edit_layer_ok():
            self.layers[index] = (material_var.get(), thickness_var.get())
            self.update_layer_table()
            dialog.destroy()
            
        ttk.Button(dialog, text='Save', command=edit_layer_ok).grid(row=2, column=0, columnspan=2, pady=10)
        
    def delete_layer(self):
        """Delete selected layer"""
        selected = self.layer_table.selection()
        if not selected:
            messagebox.showwarning("Warning", "Please select a layer to delete")
            return
            
        index = self.layer_table.index(selected[0])
        self.layers.pop(index)
        self.update_layer_table()
        
    def load_material_file(self):
        """Load material data from file"""
        name = simpledialog.askstring('Material Name', 'Enter material name:')
        if not name:
            return
            
        filepath = filedialog.askopenfilename(
            title='Select material data file',
            filetypes=[('Data files', '*.dat *.txt'), ('All files', '*.*')]
        )
        
        if filepath:
            material_files[name] = str(filepath)  # Ensure it's a string
            if name not in self.materials:
                self.materials.append(name)
            messagebox.showinfo('Success', f'Material "{name}" loaded successfully')
            
    def progress_callback(self, progress):
        """Update progress bar"""
        self.progress_var.set(progress * 100)
        self.root.update_idletasks()
        
    def calculate_3d(self):
        """Perform 3D calculation"""
        # Get parameters
        E_start = self.E_start.get()
        E_stop = self.E_stop.get()
        dE = self.dE.get()
        kx_start = self.kx_start.get()
        kx_stop = self.kx_stop.get()
        dkx = self.dkx.get()
        ky_start = self.ky_start.get()
        ky_stop = self.ky_stop.get()
        dky = self.dky.get()
        a = self.Es_amp.get()
        b = self.Ep_amp.get()
        delta_phi = self.delta_phi.get()
        
        # Estimate calculation time
        nE = int((E_stop - E_start) / dE) + 1
        nKx = int((kx_stop - kx_start) / dkx) + 1
        nKy = int((ky_stop - ky_start) / dky) + 1
        total_points = nE * nKx * nKy
        
        if total_points > 100000:
            result = messagebox.askyesno(
                "Warning", 
                f"Large calculation: {total_points:,} points.\n"
                f"This may take a very long time.\n"
                f"Consider using coarser grid for testing.\n"
                f"Continue?"
            )
            if not result:
                return
        
        self.status_label.config(text=f"Calculating {total_points:,} points...")
        self.progress_var.set(0)
        
        start_time = time.time()
        
        try:
            # Perform calculation
            self.results = calc_3d_reflectance(
                self.layers, E_start, E_stop, dE,
                kx_start, kx_stop, dkx, ky_start, ky_stop, dky,
                a, b, delta_phi, self.progress_callback
            )
            
            end_time = time.time()
            elapsed = end_time - start_time
            
            self.status_label.config(text=f"Calculation complete! Time: {elapsed:.1f}s")
            messagebox.showinfo("Success", "3D calculation completed successfully!")
            
        except Exception as e:
            self.status_label.config(text="Calculation failed")
            messagebox.showerror("Error", f"Calculation failed: {str(e)}")
            
    def visualize_results(self):
        """Visualize 3D results"""
        if self.results is None:
            messagebox.showwarning("Warning", "No results to visualize. Run calculation first.")
            return
            
        E_vals, kx_vals, ky_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = self.results
        
        # Create visualization window
        viz_window = tk.Toplevel(self.root)
        viz_window.title("3D Results Visualization")
        viz_window.geometry("1200x800")
        
        # Create notebook for different visualizations
        notebook = ttk.Notebook(viz_window)
        notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Energy slice visualization
        energy_frame = ttk.Frame(notebook)
        notebook.add(energy_frame, text="Energy Slices")
        
        # Controls for energy slice
        control_frame = ttk.Frame(energy_frame)
        control_frame.pack(fill='x', pady=5)
        
        ttk.Label(control_frame, text="Energy Index:").pack(side='left', padx=5)
        energy_var = tk.IntVar(value=len(E_vals)//2)
        energy_scale = ttk.Scale(control_frame, from_=0, to=len(E_vals)-1, 
                                variable=energy_var, orient='horizontal')
        energy_scale.pack(side='left', padx=5, fill='x', expand=True)
        
        energy_label = ttk.Label(control_frame, text=f"E = {E_vals[energy_var.get()]:.2f} eV")
        energy_label.pack(side='left', padx=5)
        
        # Quantity selection
        ttk.Label(control_frame, text="Quantity:").pack(side='left', padx=5)
        quantity_var = tk.StringVar(value="Rs")
        quantity_combo = ttk.Combobox(control_frame, textvariable=quantity_var, 
                                     values=["Rs", "Rp", "S0", "S1", "S2", "S3", "phi", "chi"], 
                                     state='readonly')
        quantity_combo.pack(side='left', padx=5)
        
        # Matplotlib canvas
        from matplotlib.figure import Figure
        fig = Figure(figsize=(10, 6))
        canvas = FigureCanvasTkAgg(fig, energy_frame)
        canvas.get_tk_widget().pack(fill='both', expand=True)
        
        def update_energy_plot():
            fig.clear()
            
            E_idx = energy_var.get()
            quantity = quantity_var.get()
            
            # Select data
            data_map = {"Rs": Rs, "Rp": Rp, "S0": S0, "S1": S1, 
                       "S2": S2, "S3": S3, "phi": phi, "chi": chi}
            data = data_map[quantity][E_idx, :, :]
            
            ax = fig.add_subplot(111)
            KX, KY = np.meshgrid(kx_vals, ky_vals, indexing='ij')
            im = ax.pcolormesh(KX, KY, data, shading='auto', cmap='viridis')
            ax.set_xlabel('kx (μm⁻¹)')
            ax.set_ylabel('ky (μm⁻¹)')
            ax.set_title(f'{quantity} at E = {E_vals[E_idx]:.2f} eV')
            fig.colorbar(im, ax=ax)
            
            energy_label.config(text=f"E = {E_vals[E_idx]:.2f} eV")
            canvas.draw()
        
        energy_scale.config(command=lambda x: update_energy_plot())
        quantity_combo.bind('<<ComboboxSelected>>', lambda x: update_energy_plot())
        
        # Initial plot
        update_energy_plot()
        
        # k-slice visualization
        k_frame = ttk.Frame(notebook)
        notebook.add(k_frame, text="k-Space Slices")
        
        # Controls for k-space slice
        k_control_frame = ttk.Frame(k_frame)
        k_control_frame.pack(fill='x', pady=5)
        
        # Slice type selection
        ttk.Label(k_control_frame, text="Slice Type:").pack(side='left', padx=5)
        slice_type_var = tk.StringVar(value="kx")
        slice_type_combo = ttk.Combobox(k_control_frame, textvariable=slice_type_var, 
                                       values=["kx", "ky"], width=8, state='readonly')
        slice_type_combo.pack(side='left', padx=5)
        
        # Slice index selection
        ttk.Label(k_control_frame, text="Slice Index:").pack(side='left', padx=5)
        k_slice_var = tk.IntVar(value=len(kx_vals)//2)
        k_slice_scale = ttk.Scale(k_control_frame, from_=0, to=len(kx_vals)-1, 
                                 variable=k_slice_var, orient='horizontal')
        k_slice_scale.pack(side='left', padx=5, fill='x', expand=True)
        
        k_slice_label = ttk.Label(k_control_frame, text="")
        k_slice_label.pack(side='left', padx=5)
        
        # Quantity selection for k-space
        ttk.Label(k_control_frame, text="Quantity:").pack(side='left', padx=5)
        k_quantity_var = tk.StringVar(value="Rs")
        k_quantity_combo = ttk.Combobox(k_control_frame, textvariable=k_quantity_var, 
                                       values=["Rs", "Rp", "S0", "S1", "S2", "S3", "phi", "chi"], 
                                       state='readonly')
        k_quantity_combo.pack(side='left', padx=5)
        
        # Matplotlib canvas for k-space
        from matplotlib.figure import Figure
        k_fig = Figure(figsize=(10, 6))
        k_canvas = FigureCanvasTkAgg(k_fig, k_frame)
        k_canvas.get_tk_widget().pack(fill='both', expand=True)
        
        def update_k_plot():
            k_fig.clear()
            
            slice_type = slice_type_var.get()
            slice_idx = k_slice_var.get()
            quantity = k_quantity_var.get()
            
            # Select data
            data_map = {"Rs": Rs, "Rp": Rp, "S0": S0, "S1": S1, 
                       "S2": S2, "S3": S3, "phi": phi, "chi": chi}
            data = data_map[quantity]
            
            ax = k_fig.add_subplot(111)
            
            if slice_type == "kx":
                # Fixed kx, plot E vs ky
                if slice_idx < len(kx_vals):
                    plot_data = data[:, slice_idx, :]
                    KY_mesh, E_mesh = np.meshgrid(ky_vals, E_vals, indexing='ij')
                    im = ax.pcolormesh(KY_mesh.T, E_mesh.T, plot_data, shading='auto', cmap='viridis')
                    ax.set_xlabel('ky (μm⁻¹)')
                    ax.set_ylabel('Energy (eV)')
                    ax.set_title(f'{quantity} at kx = {kx_vals[slice_idx]:.2f} μm⁻¹')
                    k_slice_label.config(text=f"kx = {kx_vals[slice_idx]:.2f} μm⁻¹")
                else:
                    ax.text(0.5, 0.5, 'Invalid slice index', ha='center', va='center', transform=ax.transAxes)
            else:  # slice_type == "ky"
                # Fixed ky, plot E vs kx
                if slice_idx < len(ky_vals):
                    plot_data = data[:, :, slice_idx]
                    KX_mesh, E_mesh = np.meshgrid(kx_vals, E_vals, indexing='ij')
                    im = ax.pcolormesh(KX_mesh.T, E_mesh.T, plot_data, shading='auto', cmap='viridis')
                    ax.set_xlabel('kx (μm⁻¹)')
                    ax.set_ylabel('Energy (eV)')
                    ax.set_title(f'{quantity} at ky = {ky_vals[slice_idx]:.2f} μm⁻¹')
                    k_slice_label.config(text=f"ky = {ky_vals[slice_idx]:.2f} μm⁻¹")
                else:
                    ax.text(0.5, 0.5, 'Invalid slice index', ha='center', va='center', transform=ax.transAxes)
            
            if 'im' in locals():
                k_fig.colorbar(im, ax=ax)
            
            k_canvas.draw()
        
        def update_slice_range():
            # Update the slice scale range when slice type changes
            slice_type = slice_type_var.get()
            if slice_type == "kx":
                k_slice_scale.config(to=len(kx_vals)-1)
                k_slice_var.set(len(kx_vals)//2)
            else:  # ky
                k_slice_scale.config(to=len(ky_vals)-1)
                k_slice_var.set(len(ky_vals)//2)
            update_k_plot()
        
        # Bind events
        k_slice_scale.config(command=lambda x: update_k_plot())
        k_quantity_combo.bind('<<ComboboxSelected>>', lambda x: update_k_plot())
        slice_type_combo.bind('<<ComboboxSelected>>', lambda x: update_slice_range())
        
        # Initial plot
        update_k_plot()
        
    def save_results(self):
        """Save 3D results to files"""
        if self.results is None:
            messagebox.showwarning("Warning", "No results to save. Run calculation first.")
            return
            
        folder = filedialog.askdirectory(title="Select folder to save results")
        if not folder:
            return
            
        E_vals, kx_vals, ky_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = self.results
        
        try:
            # Save coordinate arrays
            np.savetxt(os.path.join(folder, 'E_vals.dat'), E_vals)
            np.savetxt(os.path.join(folder, 'kx_vals.dat'), kx_vals)
            np.savetxt(os.path.join(folder, 'ky_vals.dat'), ky_vals)
            
            # Save 3D data arrays (flattened with metadata)
            data_arrays = {'Rs': Rs, 'Rp': Rp, 'S0': S0, 'S1': S1, 
                          'S2': S2, 'S3': S3, 'phi': phi, 'chi': chi}
            
            for name, data in data_arrays.items():
                # Save as binary numpy file for efficiency
                np.save(os.path.join(folder, f'{name}.npy'), data)
                
            # Save metadata
            metadata = {
                'shape': Rs.shape,
                'E_range': (E_vals[0], E_vals[-1]),
                'kx_range': (kx_vals[0], kx_vals[-1]),
                'ky_range': (ky_vals[0], ky_vals[-1]),
                'layers': self.layers
            }
            
            with open(os.path.join(folder, 'metadata.txt'), 'w') as f:
                for key, value in metadata.items():
                    f.write(f"{key}: {value}\n")
                    
            messagebox.showinfo("Success", f"Results saved to {folder}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save results: {str(e)}")
            
    def run(self):
        """Start the GUI"""
        self.root.mainloop()

if __name__ == '__main__':
    app = TMM3DGui()
    app.run() 