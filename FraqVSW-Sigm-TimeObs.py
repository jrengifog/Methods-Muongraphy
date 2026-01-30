import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
#Programa Calcula:
#Las fracciones + errores de Muones que pasan los detectores y el objeto
#Número de Sigmas que diferencia a fracciones para distintos materiales
#Tiempo de observación para lograr 3 Sigmas.

#Fuente para graficar
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

#from matplotlib.lines import Line2D

# -------------------- Funciones de conteo --------------------
def obtener_conteos_sin_cambioDir(L, Mat):
    file_material = file_material = f'D:/fisica pc/documentos/Carpeta-share/Flux-Mu-60min-{L}-{Mat}-Cor/'#f'Flux-Mu-60min-{L}-{Mat}-Cor/'
    N_material = np.loadtxt(file_material + 'Coinc_All.txt')
    return N_material

def obtener_conteos(L, Mat):
    base_path = 'D:/fisica pc/documentos/Carpeta-share/'#'drive/MyDrive/Colab Notebooks/'
    file_mat = f'{base_path}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/Coinc_All.txt'
    file_air = f'{base_path}Flux-Mu-60min-{L}-Air-UnDir-COR/Coinc_All.txt'
    return np.loadtxt(file_air), np.loadtxt(file_mat)

def calcular_fraccion_y_error(k, N):
    fraccion = k / N
    error = fraccion * (1 / np.sqrt(k) + 1 / np.sqrt(N))
    return fraccion, error

def FracMatvsMat(L, Mat):
    N_total, N_pasan = obtener_conteos(L, Mat)
    k1, N1 = len(N_pasan), len(N_total)
    return calcular_fraccion_y_error(k1, N1), (k1, N1)

# -------------------- Configuración de materiales --------------------
materiales = [
    {"clave": "Air", "etiqueta": "Air", "color": "red", "var": "Air", "marker": "D"},
    {"clave": "Wa", "etiqueta": "Water", "color": "blue", "var": "Wa", "marker": "^"},
    {"clave": "Con", "etiqueta": "Concrete", "color": "orange", "var": "Con", "marker": "s"},
    {"clave": "Al", "etiqueta": "Al", "color": "black", "var": "Al", "marker": "o"},
    {"clave": "Fe", "etiqueta": "Fe", "color": "green", "var": "Fe", "marker": "v"},
    {"clave": "Pb", "etiqueta": "Pb", "color": "grey", "var": "Pb", "marker": "P"},
]
anchos = [50, 100, 150, 200, 250]
resultados = {}

# -------------------- Cálculo de fracciones --------------------
fig, ax = plt.subplots(dpi=300, figsize=(7, 5))

for mat in materiales:
    frac_list, err_list, k_list, N_list = [], [], [], []
    for L in anchos:
        (f, e), (k, N) = FracMatvsMat(L, mat["var"])
        frac_list.append(f)
        err_list.append(e)
        k_list.append(k)
        N_list.append(N)

    resultados[mat["clave"]] = {
        "frac": np.array(frac_list),
        "err": np.array(err_list),
        "k": np.array(k_list),
        "N": np.array(N_list),
    }

    ax.errorbar(
        anchos, frac_list, yerr=err_list,
        fmt=mat["marker"], color=mat["color"],
        label=mat["etiqueta"], markersize=6,
        elinewidth=2, capsize=4, linestyle='none'
    )

ax.set_xlabel('Width (cm)', fontsize=16)
ax.set_ylabel('Nf/Ni', fontsize=16)
ax.tick_params(axis='both', labelsize=12)
ax.grid(True, linestyle="--", alpha=0.6)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=12)
plt.tight_layout()
plt.subplots_adjust(right=0.75)
plt.show()

# -------------------- Significancia y su error --------------------
def calcular_sigma_y_error(k1, N1, k2, N2):
    f1 = k1 / N1
    f2 = k2 / N2
    delta = abs(f1 - f2)

    # Errores individuales
    e1 = f1*np.sqrt(1 / k1 + 1 / N1)#(1 / np.sqrt(k1) + 1 / np.sqrt(N1))
    e2 = f2*np.sqrt(1 / k2 + 1 / N2)#(1 / np.sqrt(k2) + 1 / np.sqrt(N2))
    D = np.sqrt(e1**2 + e2**2)

    # Sigma (significancia)
    sigma = delta / D

    # Derivadas parciales para propagación
    df1_dk1 = (1 / N1)# - (k1 / N1**2)
    df1_dN1 = -k1 / N1**2
    df2_dk2 = (1 / N2)# - (k2 / N2**2)
    df2_dN2 = -k2 / N2**2

    dD_de1 = e1 / D
    dD_de2 = e2 / D

    de1_dk1 = (1 / np.sqrt(k1) + 1 / np.sqrt(N1)) / N1 - (0.5 * f1 / k1**1.5)
    de1_dN1 = -f1 / N1**2 + (0.5 * f1 / N1**1.5)
    de2_dk2 = (1 / np.sqrt(k2) + 1 / np.sqrt(N2)) / N2 - (0.5 * f2 / k2**1.5)
    de2_dN2 = -f2 / N2**2 + (0.5 * f2 / N2**1.5)

    dD_dk1 = dD_de1 * de1_dk1
    dD_dN1 = dD_de1 * de1_dN1
    dD_dk2 = dD_de2 * de2_dk2
    dD_dN2 = dD_de2 * de2_dN2

    dS_dk1 = (df1_dk1 * D - delta * dD_dk1) / D**2
    dS_dN1 = (df1_dN1 * D - delta * dD_dN1) / D**2
    dS_dk2 = (-df2_dk2 * D - delta * dD_dk2) / D**2
    dS_dN2 = (-df2_dN2 * D - delta * dD_dN2) / D**2

    sigma_error = np.sqrt(
        (dS_dk1 )**2 * k1 +
        (dS_dN1 )**2 * N1 +
        (dS_dk2 )**2 * k2 +
        (dS_dN2 )**2 * N2
    )

    return sigma, sigma_error
'''def ratio_uncertainty(k, N):
    """Incertidumbre de una razón de conteos de Poisson: R = k / N"""
    R = k / N
    sigma_R = R * np.sqrt(1/k + 1/N)
    return R, sigma_R

def calcular_sigma_y_error(k1, N1, k2, N2):
    # Cálculo de razones y errores
    f1, e1 = ratio_uncertainty(k1, N1)
    f2, e2 = ratio_uncertainty(k2, N2)

    delta = abs(f1 - f2)
    D = np.sqrt(e1**2 + e2**2)
    sigma = delta / D

    # Derivadas parciales de f1 y f2
    df1_dk1 = 1 / N1
    df1_dN1 = -k1 / N1**2
    df2_dk2 = 1 / N2
    df2_dN2 = -k2 / N2**2

    # Derivadas parciales de e1 y e2 con respecto a k y N
    de1_dk1 = e1 * (0.5 / k1 - 0.5 / N1)
    de1_dN1 = e1 * (-0.5 / N1 + 0.5 * k1 / N1**2)
    de2_dk2 = e2 * (0.5 / k2 - 0.5 / N2)
    de2_dN2 = e2 * (-0.5 / N2 + 0.5 * k2 / N2**2)

    # Derivadas parciales de D = sqrt(e1^2 + e2^2)
    dD_de1 = e1 / D
    dD_de2 = e2 / D

    dD_dk1 = dD_de1 * de1_dk1
    dD_dN1 = dD_de1 * de1_dN1
    dD_dk2 = dD_de2 * de2_dk2
    dD_dN2 = dD_de2 * de2_dN2

    # Derivadas parciales de sigma = delta / D
    # Nota que delta = |f1 - f2|, derivamos considerando el signo (usamos np.sign)
    sign = np.sign(f1 - f2)

    dSigma_dk1 = (sign * df1_dk1 * D - delta * dD_dk1) / D**2
    dSigma_dN1 = (sign * df1_dN1 * D - delta * dD_dN1) / D**2
    dSigma_dk2 = (-sign * df2_dk2 * D - delta * dD_dk2) / D**2
    dSigma_dN2 = (-sign * df2_dN2 * D - delta * dD_dN2) / D**2

    # Propagación de errores para sigma
    sigma_error = np.sqrt(
        dSigma_dk1**2 * k1 +
        dSigma_dN1**2 * N1 +
        dSigma_dk2**2 * k2 +
        dSigma_dN2**2 * N2
    )

    return sigma, sigma_error'''

# -------------------- Gráfico de sigmas --------------------
fig2 = plt.figure(figsize=(20, 17), dpi=100)
gs2 = fig2.add_gridspec(3, 2, hspace=0.15, wspace=0.05)
axs2 = gs2.subplots(sharey=True, sharex=True)

for ref_mat, ax in zip(materiales, axs2.flatten()):
    clave_ref = ref_mat["clave"]
    k1_list = resultados[clave_ref]["k"]
    N1_list = resultados[clave_ref]["N"]

    for mat in [m for m in materiales if m["clave"] != clave_ref]:
        clave_cmp = mat["clave"]
        k2_list = resultados[clave_cmp]["k"]
        N2_list = resultados[clave_cmp]["N"]

        sigmas, sigma_errs = [], []
        for k1, N1, k2, N2 in zip(k1_list, N1_list, k2_list, N2_list):
            s, s_err = calcular_sigma_y_error(k1, N1, k2, N2)
            sigmas.append(s)
            sigma_errs.append(s_err)

        ax.errorbar(anchos, sigmas, yerr=sigma_errs,
                    marker=mat["marker"], color=mat["color"],
                    linestyle='none', markersize=10, capsize=2,
                    label=mat["etiqueta"])
        
    ax.set_title(ref_mat["etiqueta"], fontsize=24)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.tick_params(axis='x', labelsize=20)  # ✅ ticks eje X
    ax.tick_params(axis='y', labelsize=20)  # ✅ ticks eje Y

for ax in axs2[:, 0]:
    ax.set_ylabel("σ", fontsize=24)
for ax in axs2[2, :]:
    ax.set_xlabel("Width (cm)", fontsize=20)

fig2.legend(
    handles=[Line2D([0], [0], marker=mat["marker"], color=mat["color"],
                    linestyle='none', markersize=10, label=mat["etiqueta"])
             for mat in materiales],
    loc='center right',
    bbox_to_anchor=(1.04, 0.8),
    borderaxespad=0.,
    fontsize=24,
)

plt.show()

# --------- Función para calcular tiempos y errores ----------
def calcular_tiempos_y_errores(k1s, N1s, k2s, N2s, Ndirs, kdirs):
    tiempos = []
    errores = []
    t0_vals = []

    for i in range(len(k1s)):
        sigma, sigma_err = calcular_sigma_y_error(k1s[i], N1s[i], k2s[i], N2s[i])
        t0 = 9 * kdirs[i]/Ndirs[i]
        tiempo = t0 * 9 / (sigma**2)
        err = 9 * (2 * t0 * sigma_err) / (sigma**3)
        tiempos.append(tiempo)
        errores.append(err)
        t0_vals.append(t0)

    return np.array(tiempos), np.array(errores), np.array(t0_vals)


# ---------------- Gráfico 1: Tiempo de observación -----------------
fig1 = plt.figure(figsize=(20, 17), dpi=100)
gs1 = fig1.add_gridspec(3, 2, hspace=0.15, wspace=0.05)
axs1 = gs1.subplots(sharey=True, sharex=True)

plt.yticks(fontsize=20)
#y_min_time, y_max_time = 0.001, 10
offset_step = 0  # Pequeño desplazamiento horizontal para evitar superposición

for i, (ref_mat, ax) in enumerate(zip(materiales, axs1.flatten())):
    clave_ref = ref_mat["clave"]
    ref_k = resultados[clave_ref]["k"]
    ref_N = resultados[clave_ref]["N"]

    for j, mat in enumerate(materiales):
        clave_mat = mat["clave"]
        if clave_mat == clave_ref:
            continue

        k = resultados[clave_mat]["k"]
        N = resultados[clave_mat]["N"]

        # Calcular sigma y su error con la nueva función
        sigmas = []
        sigma_errs = []
        t0_vals = []
        tiempos = []
        tiempo_errs = []

        for k1, N1, k2, N2, L in zip(ref_k, ref_N, k, N, anchos):
            sigma, sigma_err = calcular_sigma_y_error(k1, N1, k2, N2)
            sigmas.append(sigma)
            sigma_errs.append(sigma_err)

            # Obtener t0 correctamente
            N_total = obtener_conteos(L, clave_mat)[1]  # con cambio de dirección
            N_sin_dir = obtener_conteos_sin_cambioDir(L, clave_mat)
            t0 = 9 * len(N_total)/len(N_sin_dir)

            t = t0 * 9 / sigma**2
            t_err = 18 * t0 * sigma_err / sigma**3

            t0_vals.append(t0)
            tiempos.append(t)
            tiempo_errs.append(t_err)
            

        # Desplazamiento horizontal
        offset = (j - len(materiales) / 2) * offset_step
        x_offset = [x + offset for x in anchos]

        ax.errorbar(x_offset, tiempos, yerr=tiempo_errs,
                    marker=mat["marker"], color=mat["color"],
                    linestyle='none', markersize=9, capsize=3,
                    label=mat["etiqueta"], alpha=0.8)
        
    ax.set_title(f'{ref_mat["etiqueta"]}', fontsize=24)
    ax.set_yscale('log')
    #ax.set_ylim(y_min_time, y_max_time)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.tick_params(axis='x', labelsize=20)  # ✅ ticks eje X
    ax.tick_params(axis='y', labelsize=20)  # ✅ ticks eje Y
       
# Etiquetas comunes
for ax in axs1[:, 0]:
    ax.set_ylabel(r'Time to 3 $\sigma$ (hours)', fontsize=20)
for ax in axs1[2, :]:
    ax.set_xlabel("Width (cm)", fontsize=20)

# Leyenda general
fig1.legend(
    handles=[Line2D([0], [0], marker=mat["marker"], color=mat["color"],
                    linestyle='none', markersize=10, label=mat["etiqueta"])
             for mat in materiales],
    loc='center right',
    bbox_to_anchor=(1.04, 0.8),
    borderaxespad=0.,
    fontsize=24,
)
plt.show()


def calcular_tiempo_entre_Al_y_Pb(L=250):
    # Extraer datos de conteo de eventos para Al y Pb
    k_Al = resultados["Al"]["k"][anchos.index(L)]
    N_Al = resultados["Al"]["N"][anchos.index(L)]
    k_Pb = resultados["Pb"]["k"][anchos.index(L)]
    N_Pb = resultados["Pb"]["N"][anchos.index(L)]

    # Calcular sigma y su error
    sigma, sigma_err = calcular_sigma_y_error(k_Al, N_Al, k_Pb, N_Pb)

    # Calcular t0 para Pb
    N_total_Pb = obtener_conteos(250, "Air")[1]
    N_sin_dir_Pb = obtener_conteos_sin_cambioDir(L, "Air")
    t0 = 9 * len(N_total_Pb)/len(N_sin_dir_Pb)

    # Calcular tiempo y su error
    t = t0 * 9 / sigma**2
    t_err = 18 * t0 * sigma_err / sigma**3

    # Mostrar resultados
    print(f"Espesor: {L} cm")
    print(f"Sigma = {sigma:.4f} ± {sigma_err:.4f}")
    print(f"t0 = {t0:.4f} horas")
    print(f"Tiempo de observación necesario para 3σ: {t:.4f} ± {t_err:.4f} horas")

    #return t, t_err, sigma, sigma_err, t0
    
calcular_tiempo_entre_Al_y_Pb(L=250)