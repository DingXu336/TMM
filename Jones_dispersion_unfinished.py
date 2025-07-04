import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, filedialog
import os

# Physical constants
c = 3e8                    # speed of light (m/s)
mu0 = 4e-7 * np.pi         # vacuum permeability
eps0 = 1/(mu0 * c**2)      # vacuum permittivity

# Material definitions
tmaterial_files = {
    'Air': None,
    'SiO2': None
}

def get_material_eps(material, w_cm):
    if material == 'Air':
        eps = np.ones((len(w_cm), 3), dtype=complex)
    elif material == 'SiO2':
        n = 1.5
        eps_val = n**2
        eps = np.full((len(w_cm), 3), eps_val, dtype=complex)
    elif material in tmaterial_files and tmaterial_files[material] is not None:
        data = np.loadtxt(tmaterial_files[material])
        w_data = data[:, 0]
        eps_x = data[:, 1] + 1j * data[:, 2]
        eps_y = data[:, 3] + 1j * data[:, 4]
        eps_z = data[:, 5] + 1j * data[:, 6]
        eps = np.zeros((len(w_cm), 3), dtype=complex)
        eps[:, 0] = np.interp(w_cm, w_data, eps_x.real) + 1j*np.interp(w_cm, w_data, eps_x.imag)
        eps[:, 1] = np.interp(w_cm, w_data, eps_y.real) + 1j*np.interp(w_cm, w_data, eps_y.imag)
        eps[:, 2] = np.interp(w_cm, w_data, eps_z.real) + 1j*np.interp(w_cm, w_data, eps_z.imag)
    else:
        raise ValueError(f"Unknown material or file missing: {material}")
    return eps.reshape(-1, 3, 1)

def build_layers(layer_info, w_cm):
    layers = []
    for name, thickness in layer_info:
        layers.append({
            'd': thickness * 1e-9,
            'eps': get_material_eps(name, w_cm)
        })
    return layers

def local_basis_rotation(kx, ky):
    k_par = np.sqrt(kx**2 + ky**2)
    if k_par == 0:
        return np.eye(3)
    s_hat = np.array([-ky, kx, 0]) / k_par
    z_hat = np.array([0, 0, 1])
    p_hat = np.cross(z_hat, s_hat)
    return np.stack([s_hat, p_hat, z_hat], axis=1)

def rotate_eps_full(eps_tensor, kx, ky):
    R = local_basis_rotation(kx, ky)
    eps_mat = np.diagflat(eps_tensor)
    return R.T @ eps_mat @ R

def interface_jones(eps1, eps2, k_par, w):
    kz1s = np.sqrt(eps1[1,1]*(w/c)**2 - k_par**2 + 0j)
    kz2s = np.sqrt(eps2[1,1]*(w/c)**2 - k_par**2 + 0j)
    kz1p = np.sqrt(eps1[0,0]*(w/c)**2 - (eps1[0,0]/eps1[2,2])*k_par**2 + 0j)
    kz2p = np.sqrt(eps2[0,0]*(w/c)**2 - (eps2[0,0]/eps2[2,2])*k_par**2 + 0j)
    Ys1 = kz1s/(mu0*w);  Yp1 = eps0*eps1[2,2]*w/kz1p
    Ys2 = kz2s/(mu0*w);  Yp2 = eps0*eps2[2,2]*w/kz2p
    Y1 = np.diag([Ys1, Yp1]); Y2 = np.diag([Ys2, Yp2])
    return np.linalg.solve(Y2+Y1, Y2-Y1)

def gamma_jones(rj, Gp, phi):
    I = np.eye(2)
    exp2 = np.exp(2j*phi)
    num = rj + Gp*exp2
    den = I + rj@(Gp*exp2)
    return np.linalg.solve(den, num)

# 2D-kx-ky calculation at fixed energy
def calc_r_2d_single(self, E_eV, qstart, qstop, dq, layer_info, a, b, delta_phi_deg):
    w_spec = E_eV * 8065.54429
    w = 2*np.pi * c * w_spec * 100
    qx = np.arange(qstart, qstop + dq, dq) * 1e6
    qy = qx.copy()
    KX, KY = np.meshgrid(qx, qy, indexing='xy')
    Kpar = np.sqrt(KX**2 + KY**2)
    layers = build_layers(layer_info, np.array([w_spec]))
    n_ifaces = len(layers) - 1
    Gmnt = np.zeros((Kpar.shape[0],Kpar.shape[1],2,2),dtype=complex)
    for iy in range(Kpar.shape[0]):
        for ix in range(Kpar.shape[1]):
            kx, ky = KX[iy,ix], KY[iy,ix]
            kp = Kpar[iy,ix]
            e1 = rotate_eps_full(layers[-2]['eps'][0,:,0], kx, ky)
            e2 = rotate_eps_full(layers[-1]['eps'][0,:,0], kx, ky)
            G = interface_jones(e1, e2, kp, w)
            for iface in range(1, n_ifaces):
                j = len(layers) - 1 - iface
                e1 = rotate_eps_full(layers[j-1]['eps'][0,:,0], kx, ky)
                e2 = rotate_eps_full(layers[j]['eps'][0,:,0], kx, ky)
                rj = interface_jones(e1, e2, kp, w)
                kz2p = np.sqrt(e2[0,0]*(w/c)**2 - (e2[0,0]/e2[2,2])*kp**2 + 0j)
                G = gamma_jones(rj, G, kz2p*layers[j]['d'])
            Gmnt[iy,ix] = G
    δ = np.deg2rad(delta_phi_deg)
    Ein_lab = np.array([-b*np.exp(1j*δ), a, 0])
    Eref_lab = np.zeros((Kpar.shape[0], Kpar.shape[1], 3), dtype=complex)
    for iy in range(Kpar.shape[0]):
        for ix in range(Kpar.shape[1]):
            R3 = local_basis_rotation(KX[iy,ix], KY[iy,ix])
            Ein_loc = R3.T @ Ein_lab
            Eout_loc = Gmnt[iy,ix] @ Ein_loc[:2]
            Eref_lab[iy,ix] = R3 @ np.array([Eout_loc[0], Eout_loc[1], 0])
    Ex = Eref_lab[:,:,0]; Ey = Eref_lab[:,:,1]
    S0 = np.abs(Ex)**2 + np.abs(Ey)**2
    S1 = np.abs(Ex)**2 - np.abs(Ey)**2
    S2 = 2 * np.real(Ex * np.conj(Ey))
    S3 = 2 * np.imag(Ex * np.conj(Ey))
    phi = 0.5 * np.arctan2(S2, S1)
    chi = 0.5 * np.arcsin(np.clip(S3/S0, -1, 1))
    Rs = np.abs(Gmnt[:,:,0,0])**2
    Rp = np.abs(Gmnt[:,:,1,1])**2
    return KX*1e-6, KY*1e-6, Rs, Rp, S0, S1, S2, S3, phi, chi

# Energy dispersion calculation at fixed kx or ky
def calc_r_dispersion(fixed_axis, fixed_val, E_start, E_stop, dE,
                      q_start, q_stop, dq, layer_info,
                      a, b, delta_phi_deg):
    E_vals = np.arange(E_start, E_stop + dE, dE)
    q_vals = np.arange(q_start, q_stop + dq, dq)
    nE, nQ = len(E_vals), len(q_vals)
    Rs = np.zeros((nE, nQ)); Rp = np.zeros((nE, nQ))
    S0 = np.zeros((nE, nQ)); S1 = np.zeros((nE, nQ))
    S2 = np.zeros((nE, nQ)); S3 = np.zeros((nE, nQ))
    phi = np.zeros((nE, nQ)); chi = np.zeros((nE, nQ))
    for i, E in enumerate(E_vals):
        w_spec = E * 8065.54429
        w = 2*np.pi * c * w_spec * 100
        layers = build_layers(layer_info, np.array([w_spec]))
        n_ifaces = len(layers) - 1
        for j, q in enumerate(q_vals):
            kx_SI = (fixed_val*1e6) if fixed_axis=='kx' else (q*1e6)
            ky_SI = (q*1e6)           if fixed_axis=='kx' else (fixed_val*1e6)
            k_par = np.sqrt(kx_SI**2 + ky_SI**2)
            # cascade jones matrices
            e1 = rotate_eps_full(layers[-2]['eps'][0,:,0], kx_SI, ky_SI)
            e2 = rotate_eps_full(layers[-1]['eps'][0,:,0], kx_SI, ky_SI)
            G = interface_jones(e1, e2, k_par, w)
            for iface in range(1, n_ifaces):
                idx = len(layers) - 1 - iface
                e1 = rotate_eps_full(layers[idx-1]['eps'][0,:,0], kx_SI, ky_SI)
                e2 = rotate_eps_full(layers[idx]['eps'][0,:,0], kx_SI, ky_SI)
                rj = interface_jones(e1, e2, k_par, w)
                kz2p = np.sqrt(e2[0,0]*(w/c)**2 - (e2[0,0]/e2[2,2])*k_par**2 + 0j)
                G = gamma_jones(rj, G, kz2p * layers[idx]['d'])
            δ = np.deg2rad(delta_phi_deg)
            Ein_lab = np.array([-b*np.exp(1j*δ), a, 0])
            R3 = local_basis_rotation(kx_SI, ky_SI)
            Ein_loc = R3.T @ Ein_lab
            Eout_loc = G @ Ein_loc[:2]
            Eref_lab = R3 @ np.array([Eout_loc[0], Eout_loc[1], 0])
            Ex, Ey = Eref_lab[0], Eref_lab[1]
            S0[i,j] = np.abs(Ex)**2 + np.abs(Ey)**2
            S1[i,j] = np.abs(Ex)**2 - np.abs(Ey)**2
            S2[i,j] = 2 * np.real(Ex * np.conj(Ey))
            S3[i,j] = 2 * np.imag(Ex * np.conj(Ey))
            phi[i,j] = 0.5 * np.arctan2(S2[i,j], S1[i,j])
            chi[i,j] = 0.5 * np.arcsin(np.clip(S3[i,j]/S0[i,j], -1, 1))
            Rs[i,j] = np.abs(G[0,0])**2
            Rp[i,j] = np.abs(G[1,1])**2
    return E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi

class TMMGui:
    def __init__(self, root):
        self.root = root
        self.root.title("TMM Reflectivity Explorer")
        self.layers = [('Air', 1e6), ('SiO2', 500), ('Air', 1e6)]
        self.materials = list(tmaterial_files.keys())

        # Mode selection
        self.mode = tk.StringVar(value='2D k-space')
        modes = ['2D k-space', 'Energy dispersion']
        ttk.Label(root, text="Mode:").pack(side=tk.LEFT)
        ttk.Combobox(root, textvariable=self.mode, values=modes, state='readonly').pack(side=tk.LEFT, padx=4)

        # Parameters
        self.energy = tk.DoubleVar(value=1.2)
        self.Estart = tk.DoubleVar(value=1.0)
        self.Estop  = tk.DoubleVar(value=2.0)
        self.dE     = tk.DoubleVar(value=0.01)
        self.kfixed = tk.StringVar(value='kx')
        self.kval   = tk.DoubleVar(value=5.0)
        self.qstart = tk.DoubleVar(value=0)
        self.qstop  = tk.DoubleVar(value=25)
        self.qstep  = tk.DoubleVar(value=0.1)
        self.vmin   = tk.DoubleVar(value=0)
        self.vmax   = tk.DoubleVar(value=1)
        self.colormap = tk.StringVar(value='viridis')
        self.in_s_amp = tk.DoubleVar(value=1.0)
        self.in_p_amp = tk.DoubleVar(value=0.0)
        self.in_phase = tk.DoubleVar(value=0.0)

        # Layout
        param_frame = tk.Frame(root)
        param_frame.pack(pady=5)
        entries = [
            ("Energy (eV)", self.energy),
            ("E start (eV)", self.Estart),
            ("E stop (eV)", self.Estop),
            ("dE (eV)", self.dE),
            ("Fixed k-axis", self.kfixed),
            ("Fixed k-val (μm⁻¹)", self.kval),
            ("q start", self.qstart),
            ("q stop", self.qstop),
            ("dq", self.qstep),
            ("Color min", self.vmin),
            ("Color max", self.vmax),
            ("Colormap", self.colormap)
        ]
        for label, var in entries:
            frame = tk.Frame(param_frame)
            tk.Label(frame, text=label).pack(side=tk.LEFT)
            tk.Entry(frame, textvariable=var, width=6).pack(side=tk.LEFT)
            frame.pack(side=tk.LEFT, padx=4)

        jones_frame = tk.Frame(root)
        jones_frame.pack(pady=5)
        for label, var in [("|Es|", self.in_s_amp), ("|Ep|", self.in_p_amp), ("Δφ (°)", self.in_phase)]:
            tk.Label(jones_frame, text=label).pack(side=tk.LEFT)
            tk.Entry(jones_frame, textvariable=var, width=6).pack(side=tk.LEFT)

        # Layer table & controls
        self.table = ttk.Treeview(root, columns=("Material","Thickness"), show='headings')
        self.table.heading("Material", text="Material")
        self.table.heading("Thickness", text="Thickness (nm)")
        self.table.pack(pady=5)
        for m,t in self.layers:
            self.table.insert('', 'end', values=(m,t))
        btn_frame = tk.Frame(root); btn_frame.pack(pady=5)
        for txt,cmd in [("Add Layer",self.add_layer),("Edit Layer",self.edit_layer),
                        ("Delete Layer",self.delete_layer),("Load Mat File",self.load_material_file)]:
            tk.Button(btn_frame, text=txt, command=cmd).pack(side=tk.LEFT, padx=2)

        tk.Button(root, text="Compute", command=self.compute).pack(pady=5)

    def add_layer(self):
        w = tk.Toplevel(self.root)
        w.title("Add Layer")
        tk.Label(w, text="Material:").grid(row=0, column=0)
        mat = tk.StringVar(value=self.materials[0])
        ttk.Combobox(w, textvariable=mat, values=self.materials, state='readonly').grid(row=0, column=1)
        tk.Label(w, text="Thickness (nm):").grid(row=1, column=0)
        th = tk.DoubleVar()
        tk.Entry(w, textvariable=th).grid(row=1, column=1)
        def ok():
            self.layers.append((mat.get(), th.get()))
            self.table.insert('', 'end', values=(mat.get(), th.get()))
            w.destroy()
        tk.Button(w, text="Add", command=ok).grid(row=2, column=0, columnspan=2)

    def edit_layer(self):
        sel = self.table.selection()
        if not sel:
            messagebox.showinfo("Edit Layer","Select a row first.")
            return
        idx = self.table.index(sel[0])
        old_m, old_t = self.layers[idx]
        w = tk.Toplevel(self.root); w.title("Edit Layer")
        tk.Label(w, text="Material:").grid(row=0, column=0)
        mat = tk.StringVar(value=old_m)
        ttk.Combobox(w, textvariable=mat, values=self.materials, state='readonly').grid(row=0, column=1)
        tk.Label(w, text="Thickness (nm):").grid(row=1, column=0)
        th = tk.DoubleVar(value=old_t)
        tk.Entry(w, textvariable=th).grid(row=1, column=1)
        def save():
            self.layers[idx] = (mat.get(), th.get())
            self.table.item(sel[0], values=(mat.get(), th.get()))
            w.destroy()
        tk.Button(w, text="Save", command=save).grid(row=2, column=0, columnspan=2)

    def delete_layer(self):
        sel = self.table.selection()
        if not sel:
            messagebox.showinfo("Delete Layer","Select a row first.")
            return
        idx = self.table.index(sel[0])
        del self.layers[idx]
        self.table.delete(sel[0])

    def load_material_file(self):
        name = simpledialog.askstring("Material Name","Enter name:")
        if not name: return
        path = filedialog.askopenfilename(title="Select dielectric file")
        if path:
            tmaterial_files[name] = path
            if name not in self.materials:
                self.materials.append(name)
            messagebox.showinfo("Loaded", f"Material '{name}' loaded.")

    def compute(self):
        mode = self.mode.get()
        a, b, dp = self.in_s_amp.get(), self.in_p_amp.get(), self.in_phase.get()
        if mode == '2D k-space':
            E = self.energy.get()
            q0, q1, dq = self.qstart.get(), self.qstop.get(), self.qstep.get()
            X, Y, Rs, Rp, S0, S1, S2, S3, phi, chi = calc_r_2d_single(
                self, E, q0, q1, dq, self.layers, a, b, dp)
        else:
            axis, kv = self.kfixed.get(), self.kval.get()
            E0, E1, dE = self.Estart.get(), self.Estop.get(), self.dE.get()
            q0, q1, dq = self.qstart.get(), self.qstop.get(), self.qstep.get()
            E_vals, q_vals, Rs, Rp, S0, S1, S2, S3, phi, chi = calc_r_dispersion(
                axis, kv, E0, E1, dE, q0, q1, dq, self.layers, a, b, dp)
            X, Y = q_vals, E_vals
        # Plotting
        plots = [(Rs,'Rs'),(Rp,'Rp'),(S0,'S0'),(S1,'S1'),(S2,'S2'),(S3,'S3'),(phi,'phi'),(chi,'chi')]
        fig, axes = plt.subplots(2,4,figsize=(18,8))
        for idx, (data, title) in enumerate(plots):
            r, c = (0,idx) if idx<4 else (1,idx-4)
            im = axes[r,c].pcolormesh(X, Y, data, shading='auto', cmap=self.colormap.get(), vmin=self.vmin.get(), vmax=self.vmax.get())
            axes[r,c].set_title(title)
            axes[r,c].set_xlabel('kx (μm⁻¹)' if mode=='2D k-space' else 'q (μm⁻¹)')
            axes[r,c].set_ylabel('ky (μm⁻¹)' if mode=='2D k-space' else 'Energy (eV)')
            fig.colorbar(im, ax=axes[r,c])
        plt.tight_layout()
        plt.show()

if __name__ == '__main__':
    root = tk.Tk()
    app = TMMGui(root)
    root.mainloop()