   
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, filedialog

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

def calc_r_2d_single(E_eV, qstart, qstop, dq, layer_info, a=1.0, b=1.0, delta_phi_deg=0.0):
    c = 3e8
    w_spec = E_eV * 8065.54429
    w = 2*np.pi * c * w_spec * 100

    qx = np.arange(qstart, qstop + dq, dq) * 1e6
    qy = qx.copy()
    KX, KY = np.meshgrid(qx, qy, indexing='xy')
    Kpar = np.sqrt(KX**2 + KY**2)

    layers = build_layers(layer_info, np.array([w_spec]))
    n_ifaces = len(layers) - 1

    def kz_p(ex, ez, k_par):
        return np.sqrt(ex*(w/c)**2 - (ex/ez)*k_par**2 + 0j)
    def kz_s(ey, k_par):
        return np.sqrt(ey*(w/c)**2 - k_par**2 + 0j)
    def rp_p(kz1, ex1, kz2, ex2):
        return (kz1*ex2 - kz2*ex1)/(kz1*ex2 + kz2*ex1)
    def rp_s(kz1, kz2):
        return (kz1 - kz2)/(kz1 + kz2)
    def gamma(rpj, g_prev, kz2, d2):
        return (rpj + g_prev*np.exp(1j*2*kz2*d2)) / (1 + rpj*g_prev*np.exp(1j*2*kz2*d2))

    Ny, Nx = Kpar.shape
    Gs = np.zeros((Ny, Nx), dtype=complex)
    Gp = np.zeros((Ny, Nx), dtype=complex)

    for iy in range(Ny):
        for ix in range(Nx):
            kx = KX[iy, ix]
            ky = KY[iy, ix]
            k_par = Kpar[iy, ix]

            eps1_mat = rotate_eps_full(layers[-2]['eps'][0,:,0], kx, ky)
            eps2_mat = rotate_eps_full(layers[-1]['eps'][0,:,0], kx, ky)
            ex1, ey1, ez1 = eps1_mat[0,0], eps1_mat[1,1], eps1_mat[2,2]
            ex2, ey2, ez2 = eps2_mat[0,0], eps2_mat[1,1], eps2_mat[2,2]

            kz1p = kz_p(ex1, ez1, k_par)
            kz2p = kz_p(ex2, ez2, k_par)
            g_p = rp_p(kz1p, ex1, kz2p, ex2)

            kz1s = kz_s(ey1, k_par)
            kz2s = kz_s(ey2, k_par)
            g_s = rp_s(kz1s, kz2s)

            for iface in range(1, n_ifaces):
                j = len(layers) - 1 - iface
                if j < 1:
                    break
                eps1_mat = rotate_eps_full(layers[j-1]['eps'][0,:,0], kx, ky)
                eps2_mat = rotate_eps_full(layers[j  ]['eps'][0,:,0], kx, ky)
                ex1, ey1, ez1 = eps1_mat[0,0], eps1_mat[1,1], eps1_mat[2,2]
                ex2, ey2, ez2 = eps2_mat[0,0], eps2_mat[1,1], eps2_mat[2,2]

                kz2p = kz_p(ex2, ez2, k_par)
                rpj = rp_p(kz_p(ex1, ez1, k_par), ex1, kz2p, ex2)
                g_p = gamma(rpj, g_p, kz2p, layers[j]['d'])

                kz2s = kz_s(ey2, k_par)
                rpj_s = rp_s(kz_s(ey1, k_par), kz2s)
                g_s = gamma(rpj_s, g_s, kz2s, layers[j]['d'])

            Gp[iy, ix] = g_p
            Gs[iy, ix] = g_s

    delta_phi = np.deg2rad(delta_phi_deg)
    E_s = a * Gs
    E_p = b * np.exp(1j * delta_phi) * Gp

    Rs = np.abs(Gs)**2
    Rp = np.abs(Gp)**2
    S0 = np.abs(E_s)**2 + np.abs(E_p)**2
    S1 = np.abs(E_s)**2 - np.abs(E_p)**2
    S2 = 2 * np.real(E_s * np.conj(E_p))
    S3 = 2 * np.imag(E_s * np.conj(E_p))

    phi = 0.5 * np.arctan2(S2, S1)
    chi = 0.5 * np.arcsin(np.clip(S3 / S0, -1, 1))

    return (KX*1e-6, KY*1e-6, Rs, Rp, S0, S1, phi, chi)

class TMMGui:
    def __init__(self, root):
        self.root = root
        self.root.title("TMM 2D Reflectivity")

        self.layers = [('SiO2', 1e6), ('Air', 500), ('SiO2', 1e6)]
        self.materials = list(material_files.keys())

        self.energy = tk.DoubleVar(value=2.0)
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

        kx_um, ky_um, Rs, Rp, S0, S1, phi, chi = calc_r_2d_single(E, q0, q1, dq, self.layers, a, b, dp)

        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        plots = [
            (Rs,  "Rs = |rₛ|²"),
            (Rp,  "Rp = |rₚ|²"),
            (S0,  "S₀ = Rs+Rp"),
            (S1,  "S₁ = Rs−Rp"),
            (phi, "ϕ (rad)"),
            (chi, "χ (rad)")
        ]
        for idx, (data, title) in enumerate(plots):
            r = 0 if idx < 3 else 1
            c = idx if idx < 3 else (idx - 3)
            im = axes[r, c].pcolormesh(kx_um, ky_um, data, shading='auto', cmap=self.colormap.get(),
                                       vmin=self.vmin.get() if idx < 4 else None,
                                       vmax=self.vmax.get() if idx < 4 else None)
            axes[r, c].set_title(title)
            axes[r, c].set_xlabel("kx (μm⁻¹)")
            axes[r, c].set_ylabel("ky (μm⁻¹)")
            fig.colorbar(im, ax=axes[r, c])

        axes[1, 0].axis('off')
        plt.tight_layout()
        plt.show()


    def visualize_selected_layer_eps(self):
        """
        Plot dielectric tensor components ε_xx and ε_yy before and after rotation for a selected layer.
        """
        idx = simpledialog.askinteger("Layer Index", f"Enter layer index (0 to {len(self.layers)-1}):")
        if idx is None or not (0 <= idx < len(self.layers)):
            return

        # Get unrotated tensor values at given energy
        E = self.energy.get()
        w_spec = E * 8065.54429
        eps_raw = get_material_eps(self.layers[idx][0], np.array([w_spec]))[0, :, 0].real
        eps_xx_un = eps_raw[0]
        eps_yy_un = eps_raw[1]

        # Build k-grid
        q_vals = np.linspace(self.qstart.get(), self.qstop.get(), 201)
        KX, KY = np.meshgrid(q_vals, q_vals, indexing='xy')
        Ny, Nx = KX.shape

        # Prepare arrays
        eps_xx_rot = np.zeros_like(KX)
        eps_yy_rot = np.zeros_like(KX)

        # Compute rotated components
        for iy in range(Ny):
            for ix in range(Nx):
                kx = KX[iy, ix]
                ky = KY[iy, ix]
                rot_eps = rotate_eps_full(eps_raw, kx, ky)
                eps_xx_rot[iy, ix] = rot_eps[0, 0].real
                eps_yy_rot[iy, ix] = rot_eps[1, 1].real

        # Plot before and after
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        titles = ["Before Rotation: ε_xx", "Before Rotation: ε_yy",
                  "After Rotation: ε_xx",  "After Rotation: ε_yy"]
        data_list = [
            np.full_like(KX, eps_xx_un),
            np.full_like(KX, eps_yy_un),
            eps_xx_rot,
            eps_yy_rot
        ]

        for ax, data, title in zip(axes.flatten(), data_list, titles):
            im = ax.pcolormesh(KX, KY, data, shading='auto', cmap='viridis')
            ax.set_title(title)
            ax.set_xlabel("kx (μm⁻¹)")
            ax.set_ylabel("ky (μm⁻¹)")
            fig.colorbar(im, ax=ax)

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
    root.mainloop()   