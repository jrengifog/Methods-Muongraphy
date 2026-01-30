import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from matplotlib.lines import Line2D
base="D:/fisica pc/documentos/Carpeta-share/"


# -------------------- FUNCIONES DE ÁNGULO --------------------

def D_mu_sh(Data):
    IDs = Data[:, :2]
    pos = Data[:, 6:9]
    D_aft = np.hstack((IDs, pos))
    mask = np.r_[False, D_aft[1:, 1] != D_aft[:-1, 1]]
    return D_aft[mask]

def readFile_Count1fromTracks(filename):
    return D_mu_sh(np.loadtxt(filename))

def match_vectors(D1, D2):
    D1_dict = {d[1]: d for d in D1}
    D2_dict = {d[1]: d for d in D2}
    common_keys = D1_dict.keys() & D2_dict.keys()
    return np.array([D2_dict[k] for k in common_keys])

def ElementPass(F1, F2, F3, F4):
    D1 = readFile_Count1fromTracks(F1)
    D2 = readFile_Count1fromTracks(F2)
    D3 = readFile_Count1fromTracks(F3)
    D4 = readFile_Count1fromTracks(F4)
    Dbef = match_vectors(D1, D2)
    Daft = match_vectors(D3, D4)
    common_ids = set(Dbef[:, 1]) & set(Daft[:, 1])
    D4b_ = np.array([d for d in Dbef if d[1] in common_ids])
    D4a_ = np.array([d for d in Daft if d[1] in common_ids])
    return D4b_, D4a_

def angle_abs(v1, v2):
    def unit(v): return v / np.linalg.norm(v) if np.linalg.norm(v) != 0 else v
    v1_u, v2_u = unit(v1), unit(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), 0, 1.0))

def AnglefromComp(D4b, D4a):
    return [angle_abs(D4b[i, 2:], D4a[i, 2:]) for i in range(len(D4b))]

def get_angles(L, Mat):
    #base = 'drive/MyDrive/Colab Notebooks/'
    files = [f'{base}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/Track-D{i}-DMmu.txt' for i in range(1, 5)]
    DAfter, DBefore = ElementPass(*files)
    angles = np.abs(np.degrees(AnglefromComp(DAfter, DBefore)))
    return angles

def gaussiana(x, A, sigma):
    return A * np.exp(-x**2 / (2 * sigma**2))

def ajustar_angulo(angles):
    bins = np.arange(0, min(15, max(5, np.ceil(np.max(angles)))) + 1, 1)
    hist, bin_edges = np.histogram(angles, bins=bins, density=True)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    try:
        popt, pcov = curve_fit(gaussiana, centers, hist, p0=[np.max(hist), np.std(angles)],
                               bounds=([0, 0.1], [np.inf, 10]), maxfev=3000)
        sigma = popt[1]
        error = np.sqrt(np.abs(pcov[1, 1]) + 0.1**2)
    except:
        sigma = np.std(angles)
        error = np.sqrt((sigma / np.sqrt(2 * len(angles)))**2 + 0.1**2)
    return sigma, error

def calcular_significancia_con_error(ref_sigma, ref_err, sigma, err):
    delta = np.abs(ref_sigma - sigma)
    err_comb = np.sqrt(ref_err**2 + err**2)
    sig = delta / err_comb
    err_sig = sig * np.sqrt((ref_err/delta)**2 + (err/delta)**2)
    return sig, err_sig

# -------------------- FUNCIONES DE FRACCIÓN --------------------

def obtener_conteos(L, Mat):
    #base = 'drive/MyDrive/Colab Notebooks/'
    file_mat = f'{base}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/Coinc_All.txt'
    file_air = f'{base}Flux-Mu-60min-{L}-Air-UnDir-COR/Coinc_All.txt'
    return np.loadtxt(file_air), np.loadtxt(file_mat)

def calcular_sigma_y_error(k1, N1, k2, N2):
    f1, f2 = k1/N1, k2/N2
    delta = abs(f1 - f2)
    #e1 = f1 * (1/np.sqrt(k1) + 1/np.sqrt(N1))
    #e2 = f2 * (1/np.sqrt(k2) + 1/np.sqrt(N2))
    e1 = f1*np.sqrt(1 / k1 + 1 / N1)#(1 / np.sqrt(k1) + 1 / np.sqrt(N1))
    e2 = f2*np.sqrt(1 / k2 + 1 / N2)#(1 / np.sqrt(k2) + 1 / np.sqrt(N2))
    
    D = np.sqrt(e1**2 + e2**2)
    sigma = delta / D

    df1_dk1, df1_dN1 = 1/N1, -k1/N1**2
    df2_dk2, df2_dN2 = 1/N2, -k2/N2**2
    dD_de1, dD_de2 = e1/D, e2/D
    de1_dk1 = (1/np.sqrt(k1) + 1/np.sqrt(N1)) / N1 - 0.5*f1 / k1**1.5
    de1_dN1 = -f1/N1**2 + 0.5*f1 / N1**1.5
    de2_dk2 = (1/np.sqrt(k2) + 1/np.sqrt(N2)) / N2 - 0.5*f2 / k2**1.5
    de2_dN2 = -f2/N2**2 + 0.5*f2 / N2**1.5
    dD_dk1 = dD_de1 * de1_dk1
    dD_dN1 = dD_de1 * de1_dN1
    dD_dk2 = dD_de2 * de2_dk2
    dD_dN2 = dD_de2 * de2_dN2
    dS_dk1 = (df1_dk1*D - delta*dD_dk1) / D**2
    dS_dN1 = (df1_dN1*D - delta*dD_dN1) / D**2
    dS_dk2 = (-df2_dk2*D - delta*dD_dk2) / D**2
    dS_dN2 = (-df2_dN2*D - delta*dD_dN2) / D**2

    sigma_err = np.sqrt(
        (dS_dk1)**2 * k1 + (dS_dN1)**2 * N1 +
        (dS_dk2)**2 * k2 + (dS_dN2)**2 * N2
    )
    return sigma, sigma_err

# -------------------- GRAFICAR Pb vs Fe --------------------

anchos = [50, 100, 150, 200, 250]
ref_clave = "Pb"
cmp_clave = "Fe"

# Definir estilo
materiales = {
    "Pb": {"color": "grey", "marker": "P", "etiqueta": "Pb"},
    "Fe": {"color": "green", "marker": "v", "etiqueta": "Fe"}
}

sigmas_ang, errs_ang = [], []
sigmas_frac, errs_frac = [], []

for L in anchos:
    # Ángulos
    s_ref, e_ref = ajustar_angulo(get_angles(L, ref_clave))
    s_cmp, e_cmp = ajustar_angulo(get_angles(L, cmp_clave))
    s1, e1 = calcular_significancia_con_error(s_ref, e_ref, s_cmp, e_cmp)
    sigmas_ang.append(s1)
    errs_ang.append(e1)

    # Fracción
    N_air_cmp, N_cmp = obtener_conteos(L, cmp_clave)
    N_air_ref, N_ref = obtener_conteos(L, ref_clave)
    k1, N1 = len(N_cmp), len(N_air_cmp)
    k2, N2 = len(N_ref), len(N_air_ref)
    s2, e2 = calcular_sigma_y_error(k1, N1, k2, N2)
    sigmas_frac.append(s2)
    errs_frac.append(e2)

# --- Gráfico ---
plt.figure(figsize=(8, 6), dpi=300)

# Pb vs Fe (ángulo)
plt.errorbar(anchos, sigmas_ang, yerr=errs_ang,
             marker=materiales["Fe"]["marker"], color=materiales["Fe"]["color"],
             linestyle='none', label="Scattering Method (Pb vs Fe)", markersize=10, capsize=4)

# Pb vs Fe (fracción)
plt.errorbar(anchos, sigmas_frac, yerr=errs_frac,
             marker=materiales["Pb"]["marker"], color=materiales["Pb"]["color"],
             linestyle='none', label="Absorption Method (Pb vs Fe)", markersize=10, capsize=4)

#plt.title("Comparación de σ entre Pb y Fe", fontsize=16)
plt.xlabel("Width (cm)", fontsize=16)
plt.ylabel("σ", fontsize=16)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()