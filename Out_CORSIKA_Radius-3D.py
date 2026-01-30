import numpy as np
import matplotlib.pyplot as plt

# Fuente
font = {
    'family': 'serif',
    'color': 'black',
    'weight': 'normal',
    'size': 16,
}

# === ARCHIVOS ===
# Reemplaza estos paths con los reales si lo vas a ejecutar
File = 'D:/fisica pc/documentos/Carpeta-share/'
DataE6 = File + 'Proton_E2-6V0.shw.bz2'

def Read_Par(Data):
    X, Y, Z = [], [], []
    CORSIKAID, shower_id = [], []
    Rad_mu, E_mu, sh_mu, R = [], [], [], []
    Nmu = 0
    Masa_Mu = 0.106  # GeV/c^2
    DRead = np.loadtxt(Data)
    
    for i in range(len(DRead)):
        Z.append(DRead[i, 6])
        CORSIKAID.append(DRead[i, 0])
        shower_id.append(DRead[i, 7])
        R.append(np.sqrt(DRead[i, 4]**2 + DRead[i, 5]**2))
        if DRead[i, 0] in [5, 6]:  # Muons
            x = DRead[i, 4] * 1e-5
            y = DRead[i, 5] * 1e-5
            X.append(x)
            Y.append(y)
            RAD3D = np.sqrt(x**2 + y**2)
            MOMT = np.linalg.norm(DRead[i, 1:4])
            Nmu += 1
            Rad_mu.append(RAD3D)
            E_mu.append(np.sqrt(MOMT**2 + Masa_Mu**2))
            sh_mu.append(DRead[i, 7])
    
    return CORSIKAID, shower_id, X, Y, Z, Nmu, Rad_mu, E_mu, sh_mu, R

# === LECTURA ===
C_ID_6, S_ID_6, X6, Y6, Z6, Nmu6, R_mu6, E_mu6, SH_mu6, RT6 = Read_Par(DataE6)

# === GRÁFICO SCATTER POSICIONES ===
plt.figure(dpi=300, figsize=(7, 5))
sc = plt.scatter(X6, Y6, c=R_mu6, cmap='gnuplot_r', s=12, edgecolors='none')
clb = plt.colorbar(sc)
clb.ax.tick_params(labelsize=16)
clb.ax.set_title('Radius (Km)', fontsize=16)
plt.xlabel("X (km)", fontdict=font)
plt.ylabel("Y (km)", fontdict=font)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()
print(Nmu6)
print(len(S_ID_6))
print(len(C_ID_6))

# === HISTOGRAMA CORSIKA ID ===
fig, ax1 = plt.subplots(dpi=300, figsize=(7, 5))
squad = [' ', '$\gamma$', '$e^+$', '$e^-$', ' ', '$\mu^+$', '$\mu^-$',
         '$\pi^0$', '$\pi^+$', '$\pi^-$', '$K^0$', '$K^+$', '$K^-$',
         '$n^0$', '$p^+$', '$p^-$']
X1 = np.arange(16)

kwargs1 = dict(histtype='step', alpha=0.9, density=False, bins=range(1, 17),
               linewidth=2, color='darkred')

ax1.hist(C_ID_6, **kwargs1)
ax1.set_xlabel("Particle type", fontdict=font)
ax1.set_ylabel("Counts", fontdict=font)
ax1.set_yscale("log")
ax1.set_xticks(X1)

# Aplicar etiquetas con desplazamiento y tamaño de fuente
ax1.set_xticklabels(squad, fontsize=14, rotation=0, ha='left')  # ha='right' desplaza hacia la izquierda
ax1.tick_params(axis='y', labelsize=14)  # Tamaño de números en eje y

ax1.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()