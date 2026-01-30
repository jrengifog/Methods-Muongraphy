import numpy as np
import matplotlib.pyplot as plt

# Fuente para graficar
font = {
    'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 16,
}

# === Cargar archivo ===
File = 'D:/fisica pc/documentos/Carpeta-share/Flux60min_DMD.shw'
Data = np.loadtxt(File)

# === Variables ===
corsika_ids = Data[:, 0].astype(int)
PX = Data[:, 1]
PY = Data[:, 2]
PZ = Data[:, 3]
MOMENTUM = np.sqrt(PY**2 + PZ**2)#(PX**2 + 

# === Filtrar muones con PZ > 0 ===
is_muon = np.isin(corsika_ids, [5, 6])
valid_pz = PZ > 0
mask = is_muon & valid_pz

# === Cálculo del ángulo cenital para muones válidos ===
zenith_angles = np.degrees(np.arccos(PZ[mask] / MOMENTUM[mask]))

# === Histograma de ángulos ===
bins = np.arange(0, 90, 1)
muon_counts, _ = np.histogram(zenith_angles, bins=bins)
bin_centers = (bins[:-1] + bins[1:]) / 2

# === Distribución teórica cos²(θ) ===
theta_rad = np.radians(bin_centers)
theory = np.cos(theta_rad)**2
theory *= muon_counts.max() / theory.max()  # Normalizar al pico de los datos

# === Graficar ===
plt.figure(figsize=(7, 5), dpi=300)
plt.plot(bin_centers, muon_counts, marker='o', linestyle='None',
         color='blue', label='Simulated Muons', markersize=4)
plt.plot(bin_centers, theory, linestyle='--', color='black',
         label=r'Theory $\propto \cos^2(\theta)$')

plt.xlabel('Zenith Angle (degrees)', fontdict=font)
plt.ylabel('Counts', fontdict=font)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()