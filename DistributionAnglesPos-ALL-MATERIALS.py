import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import norm
from scipy.optimize import curve_fit

# --- Fuente para etiquetas
font = {'family': 'serif', 'color': 'black', 'weight': 'normal', 'size': 16}

# --- Funciones de procesamiento
def D_mu_sh(Data):
    IDs = Data[:, :2]
    pos = Data[:, 6:9]
    D_aft = np.hstack((IDs, pos))
    mask = np.r_[False, D_aft[1:, 1] != D_aft[:-1, 1]]
    return D_aft[mask]

def readFile_Count1fromTracks(filename):
    return D_mu_sh(np.loadtxt(filename))

def match_vectors(D1, D2, use='D2'):
    D1_dict = {d[1]: d for d in D1}
    D2_dict = {d[1]: d for d in D2}
    common_keys = D1_dict.keys() & D2_dict.keys()
    return np.array([D1_dict[k] if use == 'D1' else D2_dict[k] for k in common_keys])

def ElementPass(F1, F2, F3, F4):
    D1, D2 = readFile_Count1fromTracks(F1), readFile_Count1fromTracks(F2)
    D3, D4 = readFile_Count1fromTracks(F3), readFile_Count1fromTracks(F4)
    Dbef, Daft = match_vectors(D1, D2), match_vectors(D3, D4)
    common_ids = set(Dbef[:, 1]) & set(Daft[:, 1])
    D4b_ = np.array([d for d in Dbef if d[1] in common_ids])
    D4a_ = np.array([d for d in Daft if d[1] in common_ids])
    return D4b_, D4a_

def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def angle_abs(v1, v2):
    return np.arccos(np.clip(np.dot(unit_vector(v1), unit_vector(v2)), 0, 1.0))

def AnglefromComp(D4b, D4a):
    return [angle_abs(D4b[i,2:], D4a[i,2:]) for i in range(len(D4b))]

# --- Lectura y ángulos
def Read_Angle(L, Mat):
    base_path = path = 'D:/fisica pc/documentos/Carpeta-share/'
    folder = f'Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    F1 = base_path + folder + 'Track-D1-DMmu.txt'
    F2 = base_path + folder + 'Track-D2-DMmu.txt'
    F3 = base_path + folder + 'Track-D3-DMmu.txt'
    F4 = base_path + folder + 'Track-D4-DMmu.txt'
    DAfter, DBefor = ElementPass(F1, F2, F3, F4)
    return np.degrees(AnglefromComp(DAfter, DBefor))

# --- Parámetros
L = 50
#Bins = np.arange(-20, 21, 1)

materiales = [
    {"clave": "Air", "etiqueta": "Air", "color": "red",    "marker": "D"},
    {"clave": "Wa",  "etiqueta": "Water", "color": "blue",   "marker": "^"},
    {"clave": "Con", "etiqueta": "Concrete", "color": "orange", "marker": "s"},
    {"clave": "Al",  "etiqueta": "Al", "color": "black", "marker": "o"},
    {"clave": "Fe",  "etiqueta": "Fe", "color": "green", "marker": "v"},
    {"clave": "Pb",  "etiqueta": "Pb", "color": "grey",  "marker": "P"},
]

# --- Lectura de datos
datos_angulos = {mat["clave"]: Read_Angle(L, mat["clave"]) for mat in materiales}

# --- Preparar gráfico
fig, ax = plt.subplots(dpi=300, figsize=(8, 6))

# --- Plot con marcas para cada material
for mat in materiales:
    data = datos_angulos[mat["clave"]]

    #-------------------
    angle_max = np.max(data)
    angle_range = min(15, max(5, np.ceil(angle_max)))
    Bins = np.arange(0, angle_range + 10, 1)
    #-------------------

    # Histogram counts NO normalizados (para Poisson)
    counts, _ = np.histogram(data, bins=Bins)

    # Histogramo normalizado (como YA lo usas)
    hist, bin_edges = np.histogram(data, bins=Bins, density=True)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
    bin_width = np.diff(bin_edges)[0]
    N = len(data)

    # --- ERROR NORMALIZADO ---
    err_norm = np.sqrt(counts) / (N * bin_width)
    # Si un bin tiene 0 counts → error 0
    err_norm = np.nan_to_num(err_norm)

    # --- Tus marcas originales ---
    ax.plot(bin_centers, hist,
            marker=mat["marker"],
            color=mat["color"],
            label=mat["etiqueta"],
            markersize=6,
            linestyle='--',
            linewidth=0.5)

    # --- NUEVO: barras de error ---
    ax.errorbar(bin_centers, hist, yerr=err_norm,
                fmt='none', ecolor=mat["color"],
                elinewidth=0.8, capsize=2, alpha=0.8)

# --- Estética
ax.set_yscale("log")
ax.set_ylim(0.000045, None)
ax.set_xlabel("Scattering Angles ($^{\circ}$)", fontdict=font)
ax.set_ylabel("Counts (normalized)", fontdict=font)
#ax.tick_params(axis='both', labelsize=14)
ax.tick_params(axis='x', labelsize=14)  # ✅ ticks eje X
ax.tick_params(axis='y', labelsize=14)  # ✅ ticks eje Y
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

# --- Leyenda
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, fontsize=12)

plt.tight_layout()
plt.show()


