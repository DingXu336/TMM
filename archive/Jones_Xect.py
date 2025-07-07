import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, filedialog
import os

c = 3e8                    # speed of light (m/s)
mu0 = 4e-7 * np.pi         # vacuum permeability
eps0 = 1/(mu0 * c**2)      # vacuum permittivity

# Fixed material dielectric constants and dynamic loading
material_files = {
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
    elif material in material_files:
        data = np.loadtxt(material_files[material])
        w_data = data[:, 0]
        eps_x = data[:, 1] + 1j * data[:, 2]
        eps_y = data[:, 3] + 1j * data[:, 4]
        eps_z = data[:, 5] + 1j * data[:, 6]
        eps = np.zeros((len(w_cm), 3), dtype=complex)
        eps[:, 0] = np.interp(w_cm, w_data, eps_x.real) + 1j*np.interp(w_cm, w_data, eps_x.imag)
        eps[:, 1] = np.interp(w_cm, w_data, eps_y.real) + 1j*np.interp(w_cm, w_data, eps_y.imag)
        eps[:, 2] = np.interp(w_cm, w_data, eps_z.real) + 1j*np.interp(w_cm, w_data, eps_z.imag)
    else:
        raise ValueError(f"Unknown material: {material}")
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

# Jones-matrix interface and cascade
def interface_jones(eps1, eps2, k_par, w):
    # normal components
    kz1s = np.sqrt(eps1[0,0]*(w/c)**2 - k_par**2 + 0j)
    kz2s = np.sqrt(eps2[0,0]*(w/c)**2 - k_par**2 + 0j)
    kz1p = np.sqrt(eps1[1,1]*(w/c)**2 - (eps1[1,1]/eps1[2,2])*k_par**2 + 0j)
    kz2p = np.sqrt(eps2[1,1]*(w/c)**2 - (eps2[1,1]/eps2[2,2])*k_par**2 + 0j)
    # admittances
    Ys1 = kz1s/(mu0*w);  Yp1 = eps0*eps1[2,2]*w/kz1p
    Ys2 = kz2s/(mu0*w);  Yp2 = eps0*eps2[2,2]*w/kz2p
    Y1 = np.diag([Ys1, Yp1]); Y2 = np.diag([Ys2, Yp2])
    # Jones R = (Y2-Y1) @ inv(Y2+Y1)
    return np.linalg.solve(Y2+Y1, Y2-Y1)

def gamma_jones(rj, G_prev, P):
    # P = diag([e^{2iφ_s}, e^{2iφ_p}])
    num = rj + P @ G_prev @ P
    den = np.eye(2) + rj @ P @ G_prev @ P
    return np.linalg.solve(den, num)

def calc_r_2d_single(self, E_eV, qstart, qstop, dq, layer_info, a, b, delta_phi_deg):
    w_spec = E_eV * 8065.54429
    w = 2*np.pi * c * w_spec * 100

    qx = np.arange(qstart, qstop + dq, dq) * 1e6
    qy = qx.copy()
    KX, KY = np.meshgrid(qx, qy, indexing='xy')
    Kpar = np.sqrt(KX**2 + KY**2)

    layers = build_layers(layer_info, np.array([w_spec]))
    n_ifaces = len(layers) - 1

    Gmnt=np.zeros((Kpar.shape[0],Kpar.shape[1],2,2),dtype=complex)

    for iy in range(Kpar.shape[0]):
        for ix in range(Kpar.shape[1]):
            kx,ky=KX[iy,ix],KY[iy,ix]
            kp=Kpar[iy,ix]

            # rotation matrix lab->local
            R3 = local_basis_rotation(kx, ky)

            # build Jones-matrix reflection by cascading interfaces
            # bottom interface
            e1=rotate_eps_full(layers[-2]['eps'][0,:,0],kx,ky)
            e2=rotate_eps_full(layers[-1]['eps'][0,:,0],kx,ky)
            G=interface_jones(e1,e2,kp,w)
            # cascade upward
            for iface in range(1,n_ifaces):
                j=len(layers)-1-iface
                e1=rotate_eps_full(layers[j-1]['eps'][0,:,0],kx,ky)
                e2=rotate_eps_full(layers[j  ]['eps'][0,:,0],kx,ky)
                rj=interface_jones(e1,e2,kp,w)
                kz_s = np.sqrt(e2[0,0]*(w/c)**2 - kp**2 + 0j)
                kz_p = np.sqrt(e2[1,1]*(w/c)**2 - (e2[1,1]/e2[2,2])*kp**2 + 0j)
                φ_s, φ_p = kz_s * layers[j]['d'], kz_p * layers[j]['d']
                phi_prop = np.diag([np.exp(1j*φ_s), np.exp(1j*φ_p)])
                G = gamma_jones(rj, G, phi_prop)
            Gmnt[iy,ix]=G

    # input field in lab frame
    δ=np.deg2rad(delta_phi_deg)
    Ein_lab = np.array([ -b*np.exp(1j*δ), a, 0 ])   # lab: x->-p, y->s, z

    # rotate to local, apply G, rotate back
    Eref_lab = np.zeros((Kpar.shape[0], Kpar.shape[1], 3), dtype=complex)
    for iy in range(Kpar.shape[0]):
        for ix in range(Kpar.shape[1]):
            R3 = local_basis_rotation(KX[iy,ix], KY[iy,ix])
            Ein_loc = R3.T @ Ein_lab
            Eout_loc2 = Gmnt[iy,ix] @ Ein_loc[:2]
            Eref_lab[iy,ix] = R3 @ np.array([Eout_loc2[0], Eout_loc2[1], 0])

    # extract Stokes
    Ex = Eref_lab[:,:,0]; Ey = Eref_lab[:,:,1]
    S0 = np.abs(Ex)**2 + np.abs(Ey)**2
    S1 = np.abs(Ex)**2 - np.abs(Ey)**2
    S2 = 2 * np.real(Ex * np.conj(Ey))
    S3 = 2 * np.imag(Ex * np.conj(Ey))
    phi = 0.5 * np.arctan2(S2, S1)
    chi = 0.5 * np.arcsin(np.clip(S3/S0, -1, 1))

    # scalar Rs,Rp from diagonal of Gmnt
    Rs = np.abs(Gmnt[:,:,0,0])**2
    Rp = np.abs(Gmnt[:,:,1,1])**2

    return (KX*1e-6, KY*1e-6, Rs, Rp, S0, S1, S2, S3, phi, chi)


class TMMGui:
    def __init__(self, root):
        self.root = root
        self.root.title("TMM 2D Reflectivity")

        self.layers = [('Air', 1e6), ('SiO2', 500), ('Air', 1e6)]
        self.materials = list(material_files.keys())

        self.energy = tk.DoubleVar(value=1.2)
        self.qstart = tk.DoubleVar(value=0)
        self.qstop  = tk.DoubleVar(value=25)
        self.qstep  = tk.DoubleVar(value=0.1)
        self.vmin   = tk.DoubleVar(value=0)
        self.vmax   = tk.DoubleVar(value=1)
        self.colormap = tk.StringVar(value='viridis')

        self.in_s_amp  = tk.DoubleVar(value=1.0)
        self.in_p_amp  = tk.DoubleVar(value=0.0)
        self.in_phase  = tk.DoubleVar(value=0.0)

        self.table = ttk.Treeview(root, columns=("Material","Thickness"), show="headings")
        self.table.heading("Material", text="Material")
        self.table.heading("Thickness", text="Thickness (nm)")
        self.table.pack()
        for m,t in self.layers:
            self.table.insert('', 'end', values=(m,t))

        param_frame = tk.Frame(root)
        param_frame.pack()
        for label, var in [
            ("Energy (eV)",    self.energy),
            ("Start q (μm⁻¹)", self.qstart),
            ("Stop q (μm⁻¹)",  self.qstop),
            ("Δq (μm⁻¹)",      self.qstep),
            ("Color min",      self.vmin),
            ("Color max",      self.vmax),
            ("Colormap",       self.colormap),
        ]:
            row = tk.Frame(param_frame)
            tk.Label(row, text=label).pack(side=tk.LEFT)
            tk.Entry(row, textvariable=var, width=8).pack(side=tk.LEFT)
            row.pack(side=tk.LEFT, padx=4)

        jones_frame = tk.Frame(root)
        jones_frame.pack(pady=5)
        for label, var in [
            ("|Es|", self.in_s_amp),
            ("|Ep|", self.in_p_amp),
            ("Δφ (°)", self.in_phase),
        ]:
            tk.Label(jones_frame, text=label).pack(side=tk.LEFT)
            tk.Entry(jones_frame, textvariable=var, width=6).pack(side=tk.LEFT)

        btn_frame = tk.Frame(root)
        btn_frame.pack(pady=5)
        tk.Button(btn_frame, text="Add Layer",        command=self.add_layer).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Edit Layer",       command=self.edit_layer).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Delete Layer",     command=self.delete_layer).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Load Material File",command=self.load_material_file).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Compute 2D",       command=self.compute).pack(side=tk.LEFT)

    def add_layer(self):
        w = tk.Toplevel(self.root)
        w.title("Add Layer")
        tk.Label(w, text="Material:").grid(row=0, column=0, sticky="e")
        mat = tk.StringVar(value=self.materials[0])
        ttk.Combobox(w, textvariable=mat, values=self.materials, state="readonly").grid(row=0, column=1)
        tk.Label(w, text="Thickness (nm):").grid(row=1, column=0, sticky="e")
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
        w = tk.Toplevel(self.root)
        w.title("Edit Layer")
        tk.Label(w, text="Material:").grid(row=0, column=0, sticky="e")
        mat = tk.StringVar(value=old_m)
        ttk.Combobox(w, textvariable=mat, values=self.materials, state="readonly").grid(row=0, column=1)
        tk.Label(w, text="Thickness (nm):").grid(row=1, column=0, sticky="e")
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
            material_files[name] = path
            if name not in self.materials:
                self.materials.append(name)
            messagebox.showinfo("Loaded", f"Material '{name}' loaded.")

    def compute(self):
        E  = self.energy.get()
        q0 = self.qstart.get()
        q1 = self.qstop.get()
        dq = self.qstep.get()
        a  = self.in_s_amp.get()
        b  = self.in_p_amp.get()
        dp = self.in_phase.get()

        kx_um, ky_um, Rs, Rp, S0, S1, S2, S3, phi, chi = calc_r_2d_single(self, E, q0, q1, dq, self.layers, a, b, dp)

        fig, axes = plt.subplots(2, 4, figsize=(18, 8))
        plots = [
            (Rs,  "Rs = |rₛ|²", self.vmin.get(), self.vmax.get()),
            (Rp,  "Rp = |rₚ|²", self.vmin.get(), self.vmax.get()),
            (S0,  "S₀", -2, 2),
            (S1,  "S1", -1, 1),
            (S2,  "S2", -1, 1),
            (S3,  "S3", -1, 1),
            (phi, "ϕ (rad)", -np.pi/2, np.pi/2),
            (chi, "χ (rad)", -np.pi/4, np.pi/4)
        ]
        self.last_results = (kx_um, ky_um, plots)

        for idx, (data, title, vmin, vmax) in enumerate(plots):
            r = 0 if idx < 4 else 1
            c = idx if idx < 4 else (idx - 4)
            im = axes[r, c].pcolormesh(kx_um, ky_um, data, shading='auto', cmap=self.colormap.get(),
                                       vmin=vmin, vmax=vmax)
            axes[r, c].set_title(title)
            axes[r, c].set_xlabel("kx (μm⁻¹)")
            axes[r, c].set_ylabel("ky (μm⁻¹)")
            fig.colorbar(im, ax=axes[r, c])

        axes[1, 0].axis('off')
        plt.tight_layout()
        plt.show()

    def visualize_selected_layer_eps(self):
        idx = simpledialog.askinteger("Layer Index", f"Enter layer index (0 to {len(self.layers)-1}):")
        if idx is None or not (0 <= idx < len(self.layers)):
            return

        E = self.energy.get()
        w_spec = E * 8065.54429
        eps_tensor = get_material_eps(self.layers[idx][0], np.array([w_spec]))[0, :, 0].real

        q_vals = np.linspace(-10, 10, 201)
        KX, KY = np.meshgrid(q_vals, q_vals, indexing='xy')
        Nx, Ny = KX.shape

        eps_xx = np.zeros_like(KX)
        eps_xy = np.zeros_like(KX)
        eps_yy = np.zeros_like(KX)
        eps_xz = np.zeros_like(KX)
        eps_zz = np.zeros_like(KX)

        for i in range(Nx):
            for j in range(Ny):
                kx = KX[i, j]
                ky = KY[i, j]
                rot_eps = rotate_eps_full(eps_tensor, kx, ky)
                eps_xx[i, j] = rot_eps[0, 0].real
                eps_xy[i, j] = rot_eps[0, 1].real
                eps_xz[i, j] = rot_eps[0, 2].real
                eps_yy[i, j] = rot_eps[1, 1].real
                eps_zz[i, j] = rot_eps[2, 2].real


        # Plot in 2x3 grid, flatten axes for iteration
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        axes_flat = axes.flatten()
        data_list = [eps_xx, eps_xy, eps_xz, eps_yy, eps_zz]
        titles = ["Re[ε_xx]", "Re[ε_xy]", "Re[ε_xz]", "Re[ε_yy]", "Re[ε_zz]"]

        for ax, data, title in zip(axes_flat, data_list, titles):
            im = ax.pcolormesh(KX, KY, data, shading='auto', cmap='viridis')
            ax.set_title(title)
            ax.set_xlabel("kx")
            ax.set_ylabel("ky")
            fig.colorbar(im, ax=ax)

        # Hide any unused subplot (e.g., the sixth)
        for ax in axes_flat[len(data_list):]:
            ax.axis('off')

        plt.tight_layout()
        plt.show()

    def save_output(self):
        if not hasattr(self, 'last_results'):
            messagebox.showinfo("Save Error", "No data to save. Run compute first.")
            return

        folder = filedialog.askdirectory(title="Select Folder to Save Outputs")
        if not folder:
            return

        kx_um, ky_um, plots = self.last_results
        for idx, (data, title, _, _) in enumerate(plots):
            fname_base = title.split("=")[0].strip().replace(" ", "_")
            np.savetxt(os.path.join(folder, f"{fname_base}.dat"), data)

            fig, ax = plt.subplots()
            im = ax.pcolormesh(kx_um, ky_um, data, shading='auto', cmap=self.colormap.get())
            ax.set_title(title)
            ax.set_xlabel("kx (μm⁻¹)")
            ax.set_ylabel("ky (μm⁻¹)")
            fig.colorbar(im, ax=ax)
            fig.tight_layout()
            fig.savefig(os.path.join(folder, f"{fname_base}.png"))
            plt.close(fig)

def visualize_rotated_eps():
    eps_tensor = np.array([2.25, 2.25, 2.25])
    q_vals = np.linspace(-10, 10, 201)
    KX, KY = np.meshgrid(q_vals, q_vals, indexing='xy')
    Nx, Ny = KX.shape

    eps_xx = np.zeros_like(KX)
    eps_xy = np.zeros_like(KX)
    eps_yy = np.zeros_like(KX)

    for i in range(Nx):
        for j in range(Ny):
            kx = KX[i, j]
            ky = KY[i, j]
            rot_eps = rotate_eps_full(eps_tensor, kx, ky)
            eps_xx[i, j] = rot_eps[0, 0].real
            eps_xy[i, j] = rot_eps[0, 1].real
            eps_yy[i, j] = rot_eps[1, 1].real

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, data, title in zip(axes, [eps_xx, eps_xy, eps_yy], ["Re[ε_xx]", "Re[ε_xy]", "Re[ε_yy"]):
        im = ax.pcolormesh(KX, KY, data, shading='auto', cmap='viridis')
        ax.set_title(title)
        ax.set_xlabel("kx")
        ax.set_ylabel("ky")
        fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    root = tk.Tk()
    app = TMMGui(root)
    tk.Button(root, text="Visualize Dielectric Tensor", command=app.visualize_selected_layer_eps).pack(pady=5)
    tk.Button(root, text="Save Output", command=app.save_output).pack(pady=5)
    root.mainloop()

