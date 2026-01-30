import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Cargar datos
File = 'D:/fisica pc/documentos/Carpeta-share/Flux60min_DMD.shw'
data = np.loadtxt(File)

# Columnas
corsika_id = data[:, 0].astype(int)
px, py, pz = data[:, 1], data[:, 2], data[:, 3]
shower_id = data[:, 7]

# Momento total
momentum = np.sqrt(px**2 + py**2 + pz**2)

# Masas
mass_dict = {
    "muon": 0.106,
    "electron": 0.000511,
    "photon": 0.0,
    "other": 0.13957
}

# Códigos CORSIKA y estilos visuales
particle_info = {
    "muon":    {"ids": [5, 6],     "color": "g", "marker": "o"},
    "electron":{"ids": [3, 4],     "color": "b", "marker": "s"},
    "photon":  {"ids": [1],        "color": "r", "marker": "D"},
    "other":   {"ids": list(range(7, 1000)), "color": "purple", "marker": "v"}
}

# Calcular energías
energies = {}
for name, info in particle_info.items():
    mask = np.isin(corsika_id, info["ids"])
    m = mass_dict[name]
    E = np.sqrt(momentum[mask]**2 + m**2) if m > 0 else momentum[mask]
    energies[name] = E

# Bins logarítmicos
bins = np.logspace(-3, 3, 100)
bin_centers = np.sqrt(bins[:-1] * bins[1:])  # media geométrica

# Graficar con marcador
plt.figure(figsize=(7, 5), dpi=300)
for name, E in energies.items():
    counts, _ = np.histogram(E, bins=bins)
    plt.plot(bin_centers, counts, marker=particle_info[name]["marker"],
             color=particle_info[name]["color"], label=name.capitalize(), markersize=4)
#linestyle='-',
# Ejes logarítmicos
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Energy (GeV)", fontsize=16)
plt.ylabel("Counts", fontsize=16)
plt.grid(True, which="both", ls="--", lw=0.3)
plt.tight_layout()
plt.legend(fontsize=12)#, title="Particles")
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

# Estadísticas
for name, E in energies.items():
    print(f"{name.capitalize()} count: {len(E)}, Mean Energy: {np.mean(E):.3f} GeV")