import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from scipy.stats import norm

def D_mu_sh(Data):
    IDs = Data[:, :2]
    pos = Data[:, 6:9]
    D_aft = np.hstack((IDs, pos))
    mask = np.r_[False, D_aft[1:, 1] != D_aft[:-1, 1]]
    return D_aft[mask]

def readFile_Count1fromTracks(filename):
    return D_mu_sh(np.loadtxt(filename))

def count_muons_D(D):
    unique, counts = np.unique(D[:, 0], return_counts=True)
    return np.column_stack((unique, counts))

def match_vectors(D1, D2, use='D2'):
    D1_dict = {d[1]: d for d in D1}
    D2_dict = {d[1]: d for d in D2}
    common_keys = D1_dict.keys() & D2_dict.keys()
    if use == 'D1':
        return np.array([D1_dict[k] for k in common_keys])
    return np.array([D2_dict[k] for k in common_keys])

def ElementPass(F1, F2, F3, F4):
    D1, D2 = readFile_Count1fromTracks(F1), readFile_Count1fromTracks(F2)
    D3, D4 = readFile_Count1fromTracks(F3), readFile_Count1fromTracks(F4)
    Dbef = match_vectors(D1, D2, use='D2')
    Daft = match_vectors(D3, D4, use='D2')
    common_ids = set(Dbef[:, 1]) & set(Daft[:, 1])
    D4b_ = np.array([d for d in Dbef if d[1] in common_ids])
    D4a_ = np.array([d for d in Daft if d[1] in common_ids])
    return D4b_, D4a_

def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def dot_product(v1, v2):
        return np.dot(v1, v2)

def angle_abs(v1, v2):
    """ Returns the angle in radians between given vectors"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    dot_p = np.clip(dot_product(v1_u, v2_u), 0, 1.0)  # Asegura que el valor esté en [-1, 1]
    return (np.arccos(dot_p))


def AnglefromComp(D4b, D4a):
    return [angle_abs(D4b[i,2:], D4a[i,2:]) for i in range(len(D4b))]

def gaussiana(x, A, sigma):
        return A * np.exp(-x**2 / (2 * sigma**2))

# Configuración básica
Carpet = 'D:/fisica pc/documentos/Carpeta-share/'
#Bins = np.linspace(-20, 20, 30)#### cambiar los bins
# Diccionario de materiales
materiales = [
    {"clave": "Air", "etiqueta": "Air", "color": "red", "var": "Air", "marker": "D"},
    {"clave": "Wa", "etiqueta": "Water", "color": "blue", "var": "Wa", "marker": "^"},
    {"clave": "Con", "etiqueta": "Concrete", "color": "orange", "var": "Con", "marker": "s"},
    {"clave": "Al", "etiqueta": "Al", "color": "black", "var": "Al", "marker": "o"},
    {"clave": "Fe", "etiqueta": "Fe", "color": "green", "var": "Fe", "marker": "v"},
    {"clave": "Pb", "etiqueta": "Pb", "color": "grey", "var": "Pb", "marker": "P"},
]

anchos = [50, 100, 150, 200, 250]

# --- Función común para cargar y calcular ángulos absolutos ---
def get_angles(L, Mat):
    path = 'D:/fisica pc/documentos/Carpeta-share/'#'drive/MyDrive/Colab Notebooks/'
    files = [f'{path}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/Track-D{i}-DMmu.txt' for i in range(1, 5)]
    #path = f'{Carpet}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    #files = [f'{path}Track-D{i}-DMmu.txt' for i in range(1, 5)]
    DAfter, DBefore = ElementPass(*files)
    angles = np.abs(np.degrees(AnglefromComp(DAfter, DBefore)))
    return angles

# --- Función para ajuste y obtención de histograma ---
def ajustar_angulo(angles):
    angle_max = np.max(angles)
    angle_range = min(15, max(5, np.ceil(angle_max)))
    bins = np.arange(0, angle_range + 1, 1)

    hist, bin_edges = np.histogram(angles, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    try:
        A0 = np.max(hist)
        sigma0 = np.std(angles)

        params, cov = curve_fit(
            gaussiana, bin_centers, hist,
            p0=[A0, sigma0],
            bounds=([0, 0.1], [np.inf, 10]),
            maxfev=3000
        )

        A, sigma = params
        err = np.sqrt(cov[1,1]**2+0.1**2)
        #Se le agrega medio grado al error.
    except Exception as e:
        #print(f"Ajuste fallido: {e}")
        #sigma = np.std(angles)
        #err = (sigma / np.sqrt(2 * len(angles)))**2+0.1**2
        sigma = 0.389# el valor del aire que si ajustanp.std(angles)
        err = 0.1#np.sqrt(sigma)
        params = [0, sigma]
        bins, hist, bin_edges = None, None, None

    return sigma, err, bins, hist, bin_edges, params

def mostrar_distribucion_y_ajuste(L, Mat):
    angles = get_angles(L, Mat)
    sigma, err, bins, hist, bin_edges, params = ajustar_angulo(angles)

    if bins is not None:
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        x_fit = np.linspace(0, np.max(bins), 500)
        y_fit = gaussiana(x_fit, *params)
        plt.plot(x_fit, y_fit, color='red', linewidth=2, label=f"Fit {Mat} (FWHM = 1.41 ± 0.10)")

    n, bins, _ = plt.hist(angles, bins=bins, density=True, alpha=0.0)#, histtype='stepfilled', color='skyblue', edgecolor='black', label='Datos')
    
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    plt.scatter(bin_centers, n, color='blue', s=20, label=f"{Mat}")
    
    #plt.title(f'{Mat}, {L} cm  σ = {sigma:.3f} ± {err:.3f}')
    plt.xlabel('Angles (°)', fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()

#mostrar_distribucion_y_ajuste(250, 'Pb')
#exit()
# Procesamiento de materiales
def process_material(L, Mat, return_data=False):
    angles = get_angles(L, Mat)
    sigma, err, bins, hist, bin_edges, params = ajustar_angulo(angles)

    #if Mat.lower() in ['air', 'wa', 'con', 'al']:
    #    err = 0.5

    if return_data:
        return sigma, err, angles, hist, bin_edges, params
    return sigma, err


# --- Gráfico resumen: σ vs espesor ---
plt.figure(dpi=300, figsize=(7, 5))
for mat in materiales:
    clave = mat["clave"]
    sigmas, errs = [], []
    for L in anchos:
        sigma, err = process_material(L, clave)
        sigmas.append(sigma)
        errs.append(err)

    plt.errorbar(
        anchos, sigmas, yerr=errs,
        fmt=mat["marker"],
        color=mat["color"],
        markersize=6,
        linewidth=2,
        capsize=4,
        capthick=2,
        elinewidth=2,
        label=mat["etiqueta"]
    )

plt.xlabel('Width (cm)', fontsize=16)
plt.ylabel('FWHM (°)', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
plt.tight_layout()
plt.show()

# Función auxiliar para calcular significancia y su error
def calcular_significancia_con_error(ref_sigma, ref_err, sigma, err):
    diff = np.abs(ref_sigma - sigma)
    err_comb = np.sqrt(ref_err**2 + err**2)
    sig = diff / err_comb
    err_sig = sig * np.sqrt((ref_err/diff)**2 + (err/diff)**2)
    return sig, err_sig


fig2 = plt.figure(figsize=(20, 17), dpi=100)
gs2 = fig2.add_gridspec(3, 2, hspace=0.15, wspace=0.05)
axs2 = gs2.subplots(sharex=True, sharey=True)

for ref_mat, ax in zip(materiales, axs2.flatten()):
    clave_ref = ref_mat["clave"]

    for mat in [m for m in materiales if m["clave"] != clave_ref]:
        clave_cmp = mat["clave"]
        sigmaT, sigma_errs = [], []

        for L in anchos:
            sigma_ref, err_ref = process_material(L, clave_ref)
            sigma_cmp, err_cmp = process_material(L, clave_cmp)

            sig, sig_err = calcular_significancia_con_error(sigma_ref, err_ref, sigma_cmp, err_cmp)
            sigmaT.append(sig)
            sigma_errs.append(sig_err)

        ax.errorbar(
            anchos, sigmaT, yerr=sigma_errs,
            marker=mat["marker"],
            color=mat["color"],
            linestyle='none',
            markersize=10,
            capsize=4,
            linewidth=2,
            elinewidth=2,
            label=mat["etiqueta"]
        )

    ax.set_title(ref_mat["etiqueta"], fontsize=24)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.tick_params(axis='x', labelsize=20)  # ✅ ticks eje X
    ax.tick_params(axis='y', labelsize=20)  # ✅ ticks eje Y

for ax in axs2[:, 0]:
    ax.set_ylabel('σ', fontsize=24)
for ax in axs2[2, :]:
    ax.set_xlabel('Width (cm)', fontsize=20)

from matplotlib.lines import Line2D
fig2.legend(
    handles=[Line2D([0], [0], marker=m["marker"], color=m["color"],
                    linestyle='none', markersize=10, label=m["etiqueta"])
             for m in materiales],
    loc='center right',
    bbox_to_anchor=(1.04, 0.8),
    borderaxespad=0.,
    fontsize=24,
)

#plt.tight_layout()
plt.show()

def obtener_conteos_sin_cambioDir(L, Mat):
    path = f'D:/fisica pc/documentos/Carpeta-share/Flux-Mu-60min-{L}-{Mat}-Cor/'
    return np.loadtxt(path + 'Coinc_All.txt')

# Obtener número de eventos totales con cambio de dirección
def obtener_conteos_con_direccion(L, Mat):
    path = f'D:/fisica pc/documentos/Carpeta-share/Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    return np.loadtxt(path + 'Coinc_All.txt')


fig2 = plt.figure(figsize=(20, 17), dpi=100)
gs2 = fig2.add_gridspec(3, 2, hspace=0.15, wspace=0.05)
axs2 = gs2.subplots(sharex=True, sharey=True)

for i, (ref_mat, ax) in enumerate(zip(materiales, axs2.flatten())):
    clave_ref = ref_mat["clave"]

    for mat in materiales:
        clave_cmp = mat["clave"]
        if clave_cmp == clave_ref:
            continue

        tiempos = []
        errores = []

        for L in anchos:
            sigma_ref, err_ref = process_material(L, clave_ref)
            sigma_cmp, err_cmp = process_material(L, clave_cmp)

            # Significancia y su error
            sig, sig_err = calcular_significancia_con_error(sigma_ref, err_ref, sigma_cmp, err_cmp)

            # t0
            N_dir = len(obtener_conteos_con_direccion(L, clave_cmp))
            N_nodir = len(obtener_conteos_sin_cambioDir(L, clave_cmp))
            t0 = 9 * N_dir / N_nodir

            # Tiempo y error
            t_obs = 9 * t0 / sig**2
            t_err = 18 * t0 * sig_err / sig**3

            tiempos.append(t_obs)
            errores.append(t_err)

        ax.errorbar(
            anchos, tiempos, yerr=errores,
            marker=mat["marker"], color=mat["color"],
            linestyle='none', markersize=9, capsize=3,
            label=mat["etiqueta"], alpha=0.8
        )

    ax.set_title(ref_mat["etiqueta"], fontsize=18)
    ax.set_yscale('log')
    ax.grid(True, linestyle=':', alpha=0.4)
    
    ax.tick_params(axis='x', labelsize=20)  # ✅ ticks eje X
    ax.tick_params(axis='y', labelsize=20)  # ✅ ticks eje Y
    
# Etiquetas comunes
for ax in axs2[:, 0]:
    ax.set_ylabel(r'Time to 3 $\sigma$ (hours)', fontsize=24)
for ax in axs2[2, :]:
    ax.set_xlabel("Width (cm)", fontsize=20)

# Leyenda general
from matplotlib.lines import Line2D
fig2.legend(
    handles=[Line2D([0], [0], marker=mat["marker"], color=mat["color"],
                    linestyle='none', markersize=10, label=mat["etiqueta"])
             for mat in materiales],
    loc='center right',
    bbox_to_anchor=(1.04, 0.8),
    borderaxespad=0.,
    fontsize=24,
)

#plt.tight_layout()
plt.show()


def calcular_tiempo_entre_materiales(L, mat1, mat2):
    # Obtener dispersión angular y errores
    sigma1, err1 = process_material(L, mat1)
    sigma2, err2 = process_material(L, mat2)

    delta_sigma = np.abs(sigma1 - sigma2)
    delta_err = np.sqrt(err1**2 + err2**2)
    sigma = delta_sigma/delta_err
    sigma_errs = sigma * np.sqrt((err_ref/delta_sigma)**2 + (err_cmp/delta_sigma)**2)#np.sqrt(sigma)

    # Obtener t0 usando el material de referencia (por ejemplo, mat2)
    N_dir = len(obtener_conteos_con_direccion(L, "Air"))
    N_nodir = len(obtener_conteos_sin_cambioDir(L, "Air"))
    t0 = 9 * N_dir / N_nodir

    # Calcular tiempo necesario y su error
    t_obs = 9 * t0 / sigma**2
    t_err = 18 * t0 * sigma_errs / sigma**3

    # Mostrar resultados
    print(f"Comparando {mat1} vs. {mat2} a {L} cm:")
    print(f"  σ_{mat1} = {sigma1:.4f} ± {err1:.4f}")
    print(f"  σ_{mat2} = {sigma2:.4f} ± {err2:.4f}")
    print(f"  |Δσ|     = {delta_sigma:.4f} ± {delta_err:.4f}")
    print(f"  t₀       = {t0:.4f} horas")
    print(f"  Tiempo necesario para 3σ: {t_obs:.4f} ± {t_err:.4f} horas")

    return t_obs, t_err, delta_sigma, delta_err, t0, N_dir, N_nodir

calcular_tiempo_entre_materiales(250, 'Pb', 'Fe')