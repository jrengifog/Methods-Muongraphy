import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import log as ln

# Gaussiana semi
def gaussiana(x, A, sigma):
    return A * np.exp(-x**2 / (2 * sigma**2))

# Ajuste de ángulos
def ajustar_angulo(angles):
    angle_max = np.max(angles)
    angle_range = min(15, max(5, np.ceil(angle_max)))
    bins = np.arange(0, angle_range + 1, 1)
    hist, bin_edges = np.histogram(angles, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    try:
        A0 = np.max(hist)
        sigma0 = np.std(angles)
        params, cov = curve_fit(gaussiana, bin_centers, hist,
                                p0=[A0, sigma0],
                                bounds=([0, 0.1], [np.inf, 10]),
                                maxfev=3000)
        A, sigma = params
        err = np.sqrt(np.abs(cov[1, 1]**2 + 0.1**2))
    except:
        sigma = np.std(angles)
        err = (sigma / np.sqrt(2 * len(angles)))**2 + 0.1**2
    return sigma, err

# Ángulo entre vectores
def unit_vector(v): return v / np.linalg.norm(v) if np.linalg.norm(v) else v
def angle_abs(v1, v2): return np.arccos(np.clip(np.dot(unit_vector(v1), unit_vector(v2)), 0, 1.0))
def AnglefromComp(D4b, D4a): return [angle_abs(D4b[i,2:], D4a[i,2:]) for i in range(len(D4b))]

# Lectura de archivos
def D_mu_sh(Data): return np.hstack((Data[:, :2], Data[:, 6:9]))[np.r_[False, Data[1:, 1] != Data[:-1, 1]]]
def readFile(filename): return D_mu_sh(np.loadtxt(filename))
def match_vectors(D1, D2): 
    D1_dict, D2_dict = {d[1]: d for d in D1}, {d[1]: d for d in D2}
    keys = D1_dict.keys() & D2_dict.keys()
    return np.array([D2_dict[k] for k in keys]), np.array([D1_dict[k] for k in keys])

def ElementPass(F1, F2, F3, F4):
    D1, D2 = readFile(F1), readFile(F2)
    D3, D4 = readFile(F3), readFile(F4)
    Dbef, Daft = match_vectors(D1, D2), match_vectors(D3, D4)
    keys = set(Dbef[0][:,1]) & set(Daft[0][:,1])
    D4b_ = np.array([d for d in Dbef[0] if d[1] in keys])
    D4a_ = np.array([d for d in Daft[0] if d[1] in keys])
    return D4b_, D4a_

# Obtener ángulos desde archivos
def get_angles(L, Mat, base_path):
    path = f'{base_path}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    files = [f'{path}Track-D{i}-DMmu.txt' for i in range(1, 5)]
    DAfter, DBefore = ElementPass(*files)
    return np.abs(np.degrees(AnglefromComp(DAfter, DBefore)))

# Función de Molière
def calculate_moliere_sigma(energy, Z, A, D, widths): 
    X0 = 716.4 * A / (D * Z * (Z + 1) * ln(287 / np.sqrt(Z)))
    theta0 = (19.2/ energy) * np.sqrt(widths / X0) * (1 + 0.038 * ln(widths / X0))
    return np.degrees(theta0)

# Plot final para Pb con modelo teórico
def graficar_pb_con_modelo(base_path):
    L_vals = [50, 100, 150, 200, 250]
    sigmas, errores = [], []
    for L in L_vals:
        angles = get_angles(L, "Pb", base_path)
        sigma, err = ajustar_angulo(angles)
        sigmas.append(sigma)
        errores.append(err)
        
    Z, A, D, E = 82, 207.2, 11.34, 17000  # Energía del muón en MeV  # Pb
    L_teo = np.linspace(50, 250, 500)
    sigma_model = calculate_moliere_sigma(E, Z, A, D, L_teo)

    plt.figure(dpi=300, figsize=(7, 5))
    plt.errorbar(L_vals, sigmas, yerr=errores, fmt='P', color='grey',
                 markersize=6, linewidth=2, capsize=4, capthick=2,
                 elinewidth=2, label='Pb (data)')
    plt.plot(L_teo, sigma_model, 'r--', linewidth=2, label='Molière Model')
    plt.xlabel('Width (cm)', fontsize=16)
    plt.ylabel('FWHM (°)', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()

# Ejecutar
graficar_pb_con_modelo("D:/fisica pc/documentos/Carpeta-share/")