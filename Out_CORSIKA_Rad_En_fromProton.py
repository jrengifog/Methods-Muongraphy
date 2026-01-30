import numpy as np
import matplotlib.pyplot as plt

# Fuente
font = {
    'family': 'serif',
    'color': 'black',
    'weight': 'normal',
    'size': 16,
}

# Rutas
File = 'D:/fisica pc/documentos/Carpeta-share/'
archivos = {
    '0.1 TeV': File + 'Proton_E2V0.shw.bz2',
    '1 TeV': File + 'Proton_E3V0.shw.bz2',
    '10 TeV': File + 'Proton_E4V0.shw.bz2',
    '0.1–100 TeV': File + 'Proton_E2-5V0.shw.bz2',
    '0.1–1000 TeV': File + 'Proton_E2-6V0.shw.bz2'
}

colores = {
    '0.1 TeV': 'blue',
    '1 TeV': 'red',
    '10 TeV': 'green',
    '0.1–100 TeV': 'purple',
    '0.1–1000 TeV': 'orange'
}

marcadores = {
    '0.1 TeV': 'o',
    '1 TeV': 's',
    '10 TeV': '^',
    '0.1–100 TeV': 'D',
    '0.1–1000 TeV': 'v'
}

# === Función para leer datos ===
def Read_Par(Data):
    X, Y, E_mu, R = [], [], [], []
    Masa_Mu = 0.106  # GeV/c^2
    DRead = np.loadtxt(Data)
    for i in range(len(DRead)):
        if DRead[i, 0] in [5, 6]:  # Muons
            x = DRead[i, 4] * 1e-5  # km
            y = DRead[i, 5] * 1e-5  # km
            RAD3D = np.sqrt(x**2 + y**2)
            MOMT = np.linalg.norm(DRead[i, 1:4])
            E = np.sqrt(MOMT**2 + Masa_Mu**2)
            X.append(x)
            Y.append(y)
            E_mu.append(E)
            R.append(RAD3D)
    return np.array(R), np.array(E_mu)

# === Leer todos los datos ===
muon_data = {}
for label, path in archivos.items():
    R, E = Read_Par(path)
    muon_data[label] = {'R': R, 'E': E}

# === Estadísticas ===
print("\nResumen de muones por energía primaria:\n")
for label, data in muon_data.items():
    N_mu = len(data['E'])
    E_avg = np.mean(data['E'])
    R_68 = np.percentile(data['R'], 68)
    radio_68_m = R_68 * 1000  # a metros
    print(f"{label}:")
    print(f"  Número de muones     : {N_mu}")
    print(f"  Energía promedio     : {E_avg:.2f} GeV")
    print(f"  Radio al 68% de muones: {R_68:.3f} km\n")
    print(f"  Radio al 68% de muones: {radio_68_m:.3f} m\n")

# === Gráfica de radios ===
fig1, ax1 = plt.subplots(dpi=300, figsize=(7, 5))
for label, data in muon_data.items():
    r_vals, r_bins = np.histogram(data['R'], bins=70)
    bin_centers = (r_bins[:-1] + r_bins[1:]) / 2
    ax1.plot(bin_centers, r_vals, linestyle=':', marker=marcadores[label],
             markersize=4, label=label, color=colores[label])
ax1.set_xlabel("Radii (km)", fontdict=font)
ax1.set_ylabel("Counts", fontdict=font)
ax1.set_yscale("log")
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend(fontsize=12)
ax1.tick_params(axis='y', labelsize=14)
ax1.tick_params(axis='x', labelsize=14)
plt.tight_layout()
plt.show()

# === Gráfica de energías ===
fig2, ax2 = plt.subplots(dpi=300, figsize=(7, 5))
for label, data in muon_data.items():
    e_vals, e_bins = np.histogram(data['E'], bins=np.logspace(-1, 4, 70))
    bin_centers = (e_bins[:-1] + e_bins[1:]) / 2
    ax2.plot(bin_centers, e_vals, linestyle=':', marker=marcadores[label],
             markersize=4, label=label, color=colores[label])
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("$E_μ$ (GeV)", fontdict=font)
ax2.set_ylabel("Counts", fontdict=font)
ax2.grid(True, linestyle='--', alpha=0.6)
ax2.legend(fontsize=12)
ax2.tick_params(axis='y', labelsize=14)
ax2.tick_params(axis='x', labelsize=14)
plt.tight_layout()
plt.show()
