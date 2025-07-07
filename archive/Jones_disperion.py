import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, filedialog
import os

# Physical constants
c = 3e8
mu0 = 4e-7 * np.pi
eps0 = 1/(mu0 * c**2)

# Material database (name -> filepath)
material_files = {'Air': None, 'SiO2': None}

# --- Material routines ---
def get_material_eps(material, w_cm):
    if material == 'Air':
        eps = np.ones((len(w_cm), 3), dtype=complex)
    elif material == 'SiO2':
        n = 1.5
        eps = np.full((len(w_cm), 3), n**2, dtype=complex)
    elif material in material_files and material_files[material]:
        data = np.loadtxt(material_files[material])
        w_data = data[:,0]
        eps = np.zeros((len(w_cm),3), dtype=complex)
        for i, comp in enumerate([1,3,5]):
            re = np.interp(w_cm, w_data, data[:,comp])
            im = np.interp(w_cm, w_data, data[:,comp+1])
            eps[:,i] = re + 1j*im
    else:
        raise ValueError(f"Unknown material or missing file: {material}")
    return eps.reshape(-1,3,1)

def build_layers(layer_info, w_cm):
    layers = []
    for name, thickness in layer_info:
        layers.append({'d': thickness*1e-9, 'eps': get_material_eps(name, w_cm)})
    return layers

# --- Jones calculus routines ---
def local_basis_rotation(kx, ky):
    k_par = np.hypot(kx, ky)
    if k_par == 0:
        return np.eye(3)
    s_hat = np.array([-ky, kx, 0]) / k_par
    z_hat = np.array([0,0,1])
    p_hat = np.cross(z_hat, s_hat)
    return np.column_stack([s_hat, p_hat, z_hat])

def rotate_eps_full(eps_vec, kx, ky):
    R = local_basis_rotation(kx, ky)
    eps_mat = np.diag(eps_vec)
    return R.T @ eps_mat @ R

def interface_jones(eps1, eps2, k_par, w):
    kz1s = np.sqrt(eps1[0,0]*(w/c)**2 - k_par**2 + 0j)
    kz2s = np.sqrt(eps2[0,0]*(w/c)**2 - k_par**2 + 0j)
    kz1p = np.sqrt(eps1[1,1]*(w/c)**2 - (eps1[1,1]/eps1[2,2])*k_par**2 + 0j)
    kz2p = np.sqrt(eps2[1,1]*(w/c)**2 - (eps2[1,1]/eps2[2,2])*k_par**2 + 0j)
    Ys1, Yp1 = kz1s/(mu0*w), eps0*eps1[2,2]*w/kz1p
    Ys2, Yp2 = kz2s/(mu0*w), eps0*eps2[2,2]*w/kz2p
    Y1, Y2 = np.diag([Ys1,Yp1]), np.diag([Ys2,Yp2])
    return np.linalg.solve(Y2+Y1, Y2-Y1)

def gamma_jones(rj, G_prev, P):
    # P is 2×2 diag([e^{2iφ_s}, e^{2iφ_p}])
    num = rj + P @ G_prev @ P
    den = np.eye(2) + rj @ P @ G_prev @ P
    return np.linalg.solve(den, num)

# --- Dispersion calculation ---
def calc_dispersion(layer_info, fixed_axis, fixed_val,
                     E_start, E_stop, q_start, q_stop,
                     a, b, delta_phi_deg):
    # Resolution (can be parameterized later)
    dE, dq = 0.01, 0.1
    E_vals = np.arange(E_start, E_stop + dE, dE)
    q_vals = np.arange(q_start, q_stop + dq, dq)
    nE, nQ = len(E_vals), len(q_vals)
    # Allocate
    Rs = np.zeros((nE,nQ)); Rp = np.zeros((nE,nQ))
    S0 = np.zeros((nE,nQ)); S1 = np.zeros((nE,nQ))
    S2 = np.zeros((nE,nQ)); S3 = np.zeros((nE,nQ))
    phi = np.zeros((nE,nQ)); chi = np.zeros((nE,nQ))

    for i, E in enumerate(E_vals):
        w_spec = E * 8065.54429
        w = 2*np.pi * c * w_spec * 100
        layers = build_layers(layer_info, np.array([w_spec]))
        for j, q in enumerate(q_vals):
            kx = fixed_val*1e6 if fixed_axis=='kx' else q*1e6
            ky = q*1e6 if fixed_axis=='kx' else fixed_val*1e6
            k_par = np.hypot(kx, ky)
            # cascade two interfaces (Air->layer1->...->Air)
            G = None
            for iface in reversed(range(len(layers)-1)):
                e1 = rotate_eps_full(layers[iface]['eps'][0,:,0], kx, ky)
                e2 = rotate_eps_full(layers[iface+1]['eps'][0,:,0], kx, ky)
                rj = interface_jones(e1, e2, k_par, w)
                if iface == len(layers) - 2:
                    G = rj
                else:
                    kz_s = np.sqrt(e2[0,0]*(w/c)**2 - k_par**2 + 0j)
                    kz_p = np.sqrt(e2[1,1]*(w/c)**2 - (e2[1,1]/e2[2,2])*k_par**2 + 0j)
                    φ_s, φ_p = kz_s * layers[iface+1]['d'], kz_p * layers[iface+1]['d']
                    phi_prop = np.diag([np.exp(1j*φ_s), np.exp(1j*φ_p)])
                    G = gamma_jones(rj, G, phi_prop)
            # Input polarization
            δ = np.deg2rad(delta_phi_deg)
            Ein_lab = np.array([-b*np.exp(1j*δ), a, 0])
            R3 = local_basis_rotation(kx, ky)
            Ein_loc = R3.T @ Ein_lab
            Eout_loc = G @ Ein_loc[:2]
            Eref_lab = R3 @ np.array([Eout_loc[0], Eout_loc[1], 0])
            Ex, Ey = Eref_lab[0], Eref_lab[1]
            # Stokes
            S0[i,j] = abs(Ex)**2 + abs(Ey)**2
            S1[i,j] = (abs(Ex)**2 - abs(Ey)**2)/S0[i,j]
            S2[i,j] = 2*np.real(Ex*np.conj(Ey))/S0[i,j]
            S3[i,j] = 2*np.imag(Ex*np.conj(Ey))/S0[i,j]
            phi[i,j] = 0.5*np.arctan2(S2[i,j], S1[i,j])
            chi[i,j] = 0.5*np.arcsin(np.clip(S3[i,j]/S0[i,j], -1, 1))
            # Reflectances
            Rs[i,j] = abs(G[0,0])**2
            Rp[i,j] = abs(G[1,1])**2

    return E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi

# --- GUI for dispersion ---
class DispersionGui:
    def __init__(self):
        self.layers = [('SiO2', 1e6), ('Air', 500), ('SiO2', 1e6)]
        self.materials = list(material_files.keys())
        self.root = tk.Tk()
        self.root.title('Dispersion Calculator')

        # Parameter inputs
        frame = tk.Frame(self.root); frame.pack(padx=10, pady=10)
        tk.Label(frame, text='E start (eV)').grid(row=0, column=0)
        tk.Label(frame, text='E stop (eV)').grid(row=0, column=2)
        self.Estart = tk.DoubleVar(value=1.0)
        self.Estop  = tk.DoubleVar(value=3.0)
        tk.Entry(frame, textvariable=self.Estart, width=8).grid(row=0, column=1)
        tk.Entry(frame, textvariable=self.Estop,  width=8).grid(row=0, column=3)

        tk.Label(frame, text='q start (μm⁻¹)').grid(row=1, column=0)
        tk.Label(frame, text='q stop (μm⁻¹)').grid(row=1, column=2)
        self.qstart = tk.DoubleVar(value=0.0)
        self.qstop  = tk.DoubleVar(value=25.0)
        tk.Entry(frame, textvariable=self.qstart, width=8).grid(row=1, column=1)
        tk.Entry(frame, textvariable=self.qstop,  width=8).grid(row=1, column=3)

        tk.Label(frame, text='Fixed axis').grid(row=2, column=0)
        self.fixed_axis = tk.StringVar(value='kx')
        ttk.Combobox(frame, textvariable=self.fixed_axis, values=['kx','ky'], width=6, state='readonly').grid(row=2, column=1)
        tk.Label(frame, text='Fixed value (μm⁻¹)').grid(row=2, column=2)
        self.fixed_val = tk.DoubleVar(value=0.0)
        tk.Entry(frame, textvariable=self.fixed_val, width=8).grid(row=2, column=3)

        tk.Label(frame, text='|Es|').grid(row=3, column=0)
        self.Es = tk.DoubleVar(value=1.0)
        tk.Entry(frame, textvariable=self.Es, width=6).grid(row=3, column=1)
        tk.Label(frame, text='|Ep|').grid(row=3, column=2)
        self.Ep = tk.DoubleVar(value=0.0)
        tk.Entry(frame, textvariable=self.Ep, width=6).grid(row=3, column=3)
        tk.Label(frame, text='Δφ (°)').grid(row=3, column=4)
        self.dphi = tk.DoubleVar(value=0.0)
        tk.Entry(frame, textvariable=self.dphi, width=6).grid(row=3, column=5)

        # Layer table and controls
        self.table = ttk.Treeview(self.root, columns=('Material','Thickness'), show='headings')
        self.table.heading('Material', text='Material')
        self.table.heading('Thickness', text='Thickness (nm)')
        self.table.pack(padx=10, pady=5)
        for m,t in self.layers:
            self.table.insert('', 'end', values=(m,t))
        btns = tk.Frame(self.root); btns.pack()
        tk.Button(btns, text='Add Layer', command=self.add_layer).pack(side='left')
        tk.Button(btns, text='Edit Layer', command=self.edit_layer).pack(side='left')
        tk.Button(btns, text='Delete Layer', command=self.delete_layer).pack(side='left')
        tk.Button(btns, text='Load Material File', command=self.load_material_file).pack(side='left')

        # Compute/Save
        tk.Button(self.root, text='Compute', command=self.compute).pack(pady=5)
        tk.Button(self.root, text='Save Data', command=self.save_data).pack(pady=2)

        self.root.mainloop()

    # Layer management callbacks
    def add_layer(self):
        w = tk.Toplevel(self.root); w.title('Add Layer')
        tk.Label(w, text='Material').grid(row=0,column=0)
        mat = tk.StringVar(value=self.materials[0])
        ttk.Combobox(w, textvariable=mat, values=self.materials, state='readonly').grid(row=0,column=1)
        tk.Label(w, text='Thickness (nm)').grid(row=1,column=0)
        th = tk.DoubleVar()
        tk.Entry(w, textvariable=th).grid(row=1,column=1)
        def ok():
            self.layers.append((mat.get(), th.get()))
            self.table.insert('', 'end', values=(mat.get(), th.get()))
            w.destroy()
        tk.Button(w, text='Add', command=ok).grid(row=2,column=0,columnspan=2)

    def edit_layer(self):
        sel = self.table.selection()
        if not sel: return
        idx = self.table.index(sel[0])
        old_m, old_t = self.layers[idx]
        w = tk.Toplevel(self.root); w.title('Edit Layer')
        tk.Label(w, text='Material').grid(row=0,column=0)
        mat = tk.StringVar(value=old_m)
        ttk.Combobox(w, textvariable=mat, values=self.materials, state='readonly').grid(row=0,column=1)
        tk.Label(w, text='Thickness (nm)').grid(row=1,column=0)
        th = tk.DoubleVar(value=old_t)
        tk.Entry(w, textvariable=th).grid(row=1,column=1)
        def save():
            self.layers[idx] = (mat.get(), th.get())
            self.table.item(sel[0], values=(mat.get(),th.get()))
            w.destroy()
        tk.Button(w, text='Save', command=save).grid(row=2,column=0,columnspan=2)

    def delete_layer(self):
        sel = self.table.selection()
        if not sel: return
        idx = self.table.index(sel[0])
        self.layers.pop(idx)
        self.table.delete(sel[0])

    def load_material_file(self):
        name = simpledialog.askstring('Name','Material name:')
        if not name: return
        path = filedialog.askopenfilename(title='Select dielectric file')
        if path:
            material_files[name] = path
            if name not in self.materials:
                self.materials.append(name)
            messagebox.showinfo('Loaded', f"Material '{name}' loaded.")

    def compute(self):
        E0, E1 = self.Estart.get(), self.Estop.get()
        q0, q1 = self.qstart.get(), self.qstop.get()
        axis, val = self.fixed_axis.get(), self.fixed_val.get()
        a, b = self.Es.get(), self.Ep.get()
        dphi = self.dphi.get()
        self.results = calc_dispersion(self.layers, axis, val, E0, E1, q0, q1, a, b, dphi)
        E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = self.results
        # Plot
        fig, axes = plt.subplots(2,4, figsize=(16,8))
        titles = ['Rs','Rp','S0','S1','S2','S3','phi','chi']
        data = [Rs,Rp,S0,S1,S2,S3,phi,chi]
        vmin_vals = [0.0, 0.0, -2.0, -1.0, -1.0, -1.0, -np.pi/2, -np.pi/4]
        vmax_vals = [1.0, 1.0,  2.0,  1.0,  1.0,  1.0,  np.pi/2,  np.pi/4]
        for idx,(ax,title,d) in enumerate(zip(axes.flatten(), titles, data)):
            im = ax.pcolormesh(q_vals, E_vals, d, shading='auto', cmap='viridis', vmin=vmin_vals[idx], vmax=vmax_vals[idx])
            ax.set_title(title); ax.set_xlabel('q (µm⁻¹)'); ax.set_ylabel('E (eV)')
            fig.colorbar(im, ax=ax)
        plt.tight_layout(); plt.show()

    def save_data(self):
        if not hasattr(self, 'results'):
            messagebox.showinfo('Error','Run Compute first.')
            return
        E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = self.results
        folder = filedialog.askdirectory(title='Select save folder')
        if not folder: return
        np.savetxt(os.path.join(folder,'E_vals.dat'), E_vals)
        np.savetxt(os.path.join(folder,'q_vals.dat'), q_vals)
        for name, arr in zip(['Rs','Rp','S0','S1','S2','S3','phi','chi'], [Rs,Rp,S0,S1,S2,S3,phi,chi]):
            np.savetxt(os.path.join(folder,f'{name}.dat'), arr)
        messagebox.showinfo('Saved','Data files written.')

if __name__=='__main__':
    DispersionGui()
