import numpy as np
import matplotlib.pyplot as plt

# Carpeta raíz
Carpet = 'D:/fisica pc/documentos/Carpeta-share/'

# ---------------------- FUNCIONES ----------------------
def obtener_conteos(L, Mat, metodo="1A"):
    """
    Carga los datos según el método elegido.
    metodo = "1A" usa coincidences_Bef.txt
    metodo = "1B" usa Coinc_All.txt
    """
    file_material = f'Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    file_aire = f'Flux-Mu-60min-{L}-Air-UnDir-COR/'  # Air
    
    if metodo == "1A":
        N_aire = np.loadtxt(Carpet + file_aire + 'Coinc_All.txt')
    elif metodo == "1B":
        N_aire = np.loadtxt(Carpet + file_aire + 'coincidences_Bef.txt')
    else:
        raise ValueError("Método no reconocido. Usa '1A' o '1B'")
    
    N_material = np.loadtxt(Carpet + file_material + 'Coinc_All.txt')
    return N_aire, N_material

def calcular_fraccion(N_total, N_pasan):
    fraccion = len(N_pasan) / len(N_total)
    error = ((np.sqrt(len(N_total)) / len(N_total)) +
             (np.sqrt(len(N_pasan)) / len(N_pasan))) * fraccion
    return fraccion, error

def Fraq_por_material(Mat, anchos, metodo="1A"):
    """
    Calcula fracciones para un material y un método.
    """
    FracM, ErF = [], []
    for L in anchos:
        N_aire, N_material = obtener_conteos(L, Mat, metodo=metodo)
        Frac, er = calcular_fraccion(N_aire, N_material)
        FracM.append(Frac)
        ErF.append(er)
    return FracM, ErF

# ---------------------- CONFIG ----------------------
M_Fe = "Fe"
W_Mat = [50, 100, 150, 200, 250]

# ---------------------- CÁLCULOS ----------------------
frac_1A, err_1A = Fraq_por_material(M_Fe, W_Mat, metodo="1A")
frac_1B, err_1B = Fraq_por_material(M_Fe, W_Mat, metodo="1B")

# ---------------------- GRAFICAR ----------------------
fig, ax = plt.subplots(dpi=300, figsize=(7, 5))

ax.errorbar(W_Mat, frac_1A, yerr=err_1A,
            fmt="o", color="green", ecolor="green",
            elinewidth=2, capsize=4, markersize=6,
            label="Fe - D4/D1 Method 1A")

ax.errorbar(W_Mat, frac_1B, yerr=err_1B,
            fmt="s", color="blue", ecolor="blue",
            elinewidth=2, capsize=4, markersize=6,
            label="Fe - D4/D4(Air) Method 1B")

# ---------------------- ESTILO ----------------------
ax.set_xlabel("Width (cm)", fontsize=16)
ax.set_ylabel("Nf/Ni", fontsize=16)
ax.tick_params(axis="both", labelsize=12)
ax.grid(True, linestyle="--", alpha=0.6)
ax.legend(fontsize=12)
plt.tight_layout()
plt.show()
