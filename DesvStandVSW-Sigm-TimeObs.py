#Ayuda para mejorar codigo por partes y optimizar ChatGPT
#Programa Calcula:
#Distribucion_Angulos-Desv.Estandars + errores
#Número de Sigmas que diferencia a Distribucion_Angulos-Desv.Estandars para distintos materiales
#Tiempo de observación para lograr 3 Sigmas.
#TOMANDO EN CUENTA ANGULOS NEGATIVOS.
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def D_mu_sh(Data):
    IDs = Data[:, :2]
    pos = Data[:, 6:9]
    D_aft = np.hstack((IDs, pos))
    mask = np.r_[False, D_aft[1:, 1] != D_aft[:-1, 1]]
    return D_aft[mask]

def readFile_Count1fromTracks(filename):
    return D_mu_sh(np.loadtxt(filename))

#def count_muons_D(D):
#    unique, counts = np.unique(D[:, 0], return_counts=True)
#    return np.column_stack((unique, counts))

def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def angle_abs_deg(D4b, D4a):
    v1, v2 = D4a[:, 2:], D4b[:, 2:]
    dot = np.einsum('ij,ij->i', v1, v2)
    norms = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        cos_theta = np.clip(np.true_divide(dot, norms), -1.0, 1.0)
        cos_theta[np.isnan(cos_theta)] = 1.0
    return np.degrees(np.arccos(cos_theta))

#def AnglefromComp(D4b, D4a):
#    return angle_abs_deg(D4b, D4a).tolist()

def angle_s(v1, v2):
    v1_u, v2_u = unit_vector(v1), unit_vector(v2)
    det2D = np.linalg.det([v1_u[-2:], v2_u[-2:]])
    sign = 1 if det2D == 0 else -np.sign(det2D)
    dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return sign * np.arccos(dot)
#El signo me permite realizar el ajuste gausiano con la distribución de ángulos.

def AnglefromComp(D4b, D4a):
    return [angle_s(D4b[i,2:], D4a[i,2:]) for i in range(len(D4b))]



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

Carpet = 'D:/fisica pc/documentos/Carpeta-share/'
Bins = np.arange(-25, 26, 1)
def gaussian(x, mean, amplitude, stddev):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def Read_Angle(L, Mat):
    path = f'{Carpet}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    files = [path + f'Track-D{i}-DMmu.txt' for i in range(1, 5)]
    DAfter, DBefore = ElementPass(*files)
    angles = np.degrees(AnglefromComp(DAfter, DBefore))

    if len(angles) < 10:
        return 0.0, 1.0  # Retornar valor neutro si hay pocos datos
    #angles = angles[np.abs(angles) < 20]  # filtro anti-outlier
    #descartar ángulos extremos que probablemente sean errores numéricos o outliers (p. ej. > 20°):
        
    bin_heights, bin_borders = np.histogram(angles, Bins, density=True)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

    try:
        par, cov = curve_fit(
            gaussian, bin_centers, bin_heights,
            p0=[0., 1., 1.], maxfev=2000
        )
        err = np.sqrt(cov[2, 2])# + 0.05 ** 2)
        return par[2], err
    except Exception as e:
        print(f'Error en ajuste Gaussiano: {e}')
        return 0.0, 1.0  # Valor de reserva

# Comparación ejemplo
#SiAir_W = NSIG(*results['Al'], *results['Pb'])
#print(SiAir_W)

def mostrar_distribucion_y_ajuste(L, Mat):
    path = f'{Carpet}Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    files = [path + f'Track-D{i}-DMmu.txt' for i in range(1, 5)]
    DAfter, DBefore = ElementPass(*files)
    angles = np.degrees(AnglefromComp(DAfter, DBefore))


    if len(angles) < 10:
        print("Muy pocos ángulos para mostrar distribución.")
        return

    bin_heights, bin_borders = np.histogram(angles, Bins, density=True)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

    try:
        par, cov = curve_fit(gaussian, bin_centers, bin_heights, p0=[0., 1., 1.], maxfev=2000)
        x_fit = np.linspace(-25, 25, 500)
        y_fit = gaussian(x_fit, *par)
        plt.plot(x_fit, y_fit, color='red', linewidth=2, label='Ajuste Gaussiano')
        stddev_err = np.sqrt(cov[2, 2])
    except Exception as e:
        print(f'Error en ajuste: {e}')
        par = [0, 0, 0]
        stddev_err = 0

    plt.hist(angles, Bins, density=True, alpha=0.5, histtype='stepfilled',
             color='skyblue', edgecolor='black', label='Datos')
    plt.title(f'{Mat}, {L} cm $\sigma$ = {par[2]:.3f} ± {stddev_err:.3f}')
    plt.xlabel('Ángulo (grados)')
    plt.ylabel('Densidad')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    mean, amplitude, stddev = par
    #print(f'  Media     = {mean:.3f}')
    #print(f'  Amplitud  = {amplitude:.3f}')
    #print(f'  Std. Dev  = {stddev:.3f} ± {stddev_err:.3f}')

# Ejecutar ejemplo
mostrar_distribucion_y_ajuste(200, 'Wa')

def gaussian(x, mean, amplitude, stddev):
    """Función gaussiana estándar para ajuste."""
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))


def obtener_angulos(L, material, aplicar_filtro=True):
    """Lee archivos de trayectoria y calcula ángulos de dispersión."""
    path = f'{Carpet}Flux-Mu-60min-{L}-{material}-UnDir-COR/'
    archivos = [path + f'Track-D{i}-DMmu.txt' for i in range(1, 5)]
    DAfter, DBefore = ElementPass(*archivos)
    angulos = np.degrees(AnglefromComp(DAfter, DBefore))

    if aplicar_filtro:
        angulos = angulos[np.abs(angulos) < 20]  # filtro anti-outliers

    return angulos


def ajustar_gaussiana(angulos):
    """Ajusta una gaussiana a la distribución de ángulos."""
    if len(angulos) < 10:
        return (0.0, 1.0), None, None

    bin_heights, bin_borders = np.histogram(angulos, Bins, density=True)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

    try:
        parametros, cov = curve_fit(
            gaussian, bin_centers, bin_heights,
            p0=[0., 1., 1.], maxfev=3000
        )
        error_stddev = np.sqrt(cov[2, 2]**2+0.1**2)
        return (parametros[2], error_stddev), parametros, cov
    except Exception as e:
        print(f'Error en ajuste gaussiano: {e}')
        return (0.0, 1.0), None, None


'''def graficar_distribucion(angulos, parametros=None, L=None, material=None):
    """Grafica la distribución de ángulos y el ajuste gaussiano (si disponible)."""
    bin_heights, bin_borders = np.histogram(angulos, Bins, density=True)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

    plt.hist(angulos, Bins, density=True, alpha=0.5, histtype='stepfilled',
             color='skyblue', edgecolor='black', label='Datos')

    if parametros is not None:
        x_fit = np.linspace(Bins[0], Bins[-1], 500)
        y_fit = gaussian(x_fit, *parametros)
        plt.plot(x_fit, y_fit, color='red', linewidth=2, label='Ajuste Gaussiano')
        plt.text(-20, 0.9 * max(bin_heights), f"$\sigma$ = {parametros[2]:.2f}", fontsize=10)

    plt.title(f'{material}, {L} cm $\sigma$ = {parametros[2]:.2f}')
    plt.xlabel('Ángulo (grados)')
    plt.ylabel('Densidad')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()'''


def analizar_material(L, material, graficar=True, aplicar_filtro=True):
    """Función principal: calcula y ajusta dispersión angular para un material."""
    angulos = obtener_angulos(L, material, aplicar_filtro)

    if len(angulos) < 10:
        print("Muy pocos datos para análisis.")
        return (0.0, 1.0)

    (sigma, sigma_err), parametros, _ = ajustar_gaussiana(angulos)

    #if graficar:
    #    graficar_distribucion(angulos, parametros, L, material)

    return sigma, sigma_err

#analizar_material(250, 'Pb', graficar=True)
#analizar_material(250, 'Pb', graficar=False, aplicar_filtro=False)
#angulos = obtener_angulos(250, 'Pb')
#_, parametros, _ = ajustar_gaussiana(angulos)
#graficar_distribucion(angulos, parametros, L=250, material='Pb')
#exit()

def arr_Sig(material):
    layers = [50, 100, 150, 200, 250]#
    results = [analizar_material(L, material) for L in layers]#[Read_Angle(L, material) for L in layers]
    sigs, errs = zip(*results)
    return list(sigs), list(errs)

'''def NSigmas(x, erx, y, ery):
    return [round(abs(xi - yi) / np.sqrt(ex**2 + ey**2), 3)
            for xi, yi, ex, ey in zip(x, y, erx, ery)]

def NSIG(SIG1, ERR1, SIG2, ERR2):
    return NSigmas(SIG1, ERR1, SIG2, ERR2)'''

# Espesores considerados
espesores = [50, 100, 150, 200, 250]
# Marcas y colores personalizados para cada material
marcas = ['o', 's', 'D','^',  'v', 'P']
colores = ['black', 'orange', 'red', 'blue', 'green', 'grey']
# Lista de materiales en orden
materiales = ['Al', 'Con', 'Air', 'Wa', 'Fe', 'Pb']
results = {mat: arr_Sig(mat) for mat in materiales}

fig, ax = plt.subplots(figsize=(10, 6))
# Graficar cada material
for i, mat in enumerate(materiales):
    sigmas, errs = results[mat]
    ax.errorbar(
        espesores, sigmas, yerr=errs,
        fmt=marcas[i], color=colores[i],
        label=mat, capsize=4, markersize=6, linewidth=1.5
    )

# Personalización
ax.set_xlabel('Width (cm)', fontsize=14)
ax.set_ylabel('Gaussian fit width($^{\circ}$)', fontsize=14)
#ax.set_title('Desviación estándar del ángulo vs. Espesor', fontsize=14)
ax.tick_params(axis='both', labelsize=12)
ax.grid(True, linestyle="--", alpha=0.6)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=10)

# Ajustes finales
plt.tight_layout()
plt.subplots_adjust(right=0.75)
plt.show()

print("Número de bins:", len(Bins)-1)
print("Anchura de cada bin:", np.diff(Bins)[0])

from matplotlib.lines import Line2D

'''def graficar_comparaciones_significativas():
    # Configuración de la figura
    fig, axs = plt.subplots(3, 2, figsize=(20, 17), dpi=100)
    fig.subplots_adjust(hspace=0.3, wspace=0.2, right=0.85)
    
    # Espesores considerados
    espesores = np.array([50, 100, 150, 200, 250])
    
    # Configuración visual de materiales
    config_materiales = {
        'Al': {'marker': 'o', 'color': 'black', 'label': 'Al'},
        'Con': {'marker': 's', 'color': 'orange', 'label': 'Concrete'},
        'Air': {'marker': 'D', 'color': 'red', 'label': 'Air'},
        'Wa': {'marker': '^', 'color': 'blue', 'label': 'Water'},
        'Fe': {'marker': 'v', 'color': 'green', 'label': 'Fe'},
        'Pb': {'marker': 'P', 'color': 'grey', 'label': 'Pb'}
    }
    # Determinar el rango Y máximo de todos los datos primero
    y_max = 0
    for ref_mat in materiales:
        ref_sigmas, ref_errs = results[ref_mat]
        for mat in materiales:
            if mat == ref_mat:
                continue
            sigmas, errs = results[mat]
            for i in range(len(espesores)):
                try:
                    sig = np.abs(ref_sigmas[i] - sigmas[i]) / np.sqrt(ref_errs[i]**2 + errs[i]**2)
                    if not np.isnan(sig) and sig > y_max:
                        y_max = sig
                except:
                    continue
    y_max = min(80, y_max * 1.2)  # Límite máximo de 20 con 20% de margen
    #factor_amplificacion_error = 2
    # Función para calcular significancias
    def calcular_significancia_con_error(ref_sigma, ref_err, sigma, err):
        with np.errstate(divide='ignore', invalid='ignore'):
            diff = np.abs(ref_sigma - sigma)
            err_comb = np.sqrt(ref_err**2 + err**2)
            sig = np.divide(diff, err_comb, out=np.zeros_like(diff), where=err_comb>0)
            err_sig = sig * np.sqrt((ref_err/diff)**2 + (err_comb/(ref_sigma-sigma))**2) if (diff > 0) else 0
            return sig, err_sig# * factor_amplificacion_error  # Amplificación solo visual
        
    # Generar cada subplot
    for idx, ref_mat in enumerate(materiales):
        ax = axs[idx//2, idx%2]
        
        # Obtener datos de referencia
        ref_sigmas, ref_errs = results[ref_mat]
        
        # Graficar cada comparación
        for mat in materiales:
            if mat == ref_mat:
                continue
                
            # Obtener datos del material a comparar
            sigmas, errs = results[mat]
            
            # Calcular significancias para cada espesor
            sigmas_comp = []
            errs_comp = []
            
            for i in range(len(espesores)):
                try:
                    sig, err_sig = calcular_significancia_con_error(
                        ref_sigmas[i], ref_errs[i],
                        sigmas[i], errs[i]
                    )
                    sigmas_comp.append(sig)
                    errs_comp.append(err_sig)
                except:
                    sigmas_comp.append(np.nan)
                    errs_comp.append(np.nan)
            
            # Graficar solo markers con error bars
            ax.errorbar(espesores, sigmas_comp, yerr=errs_comp,
                       fmt=config_materiales[mat]['marker'],
                       color=config_materiales[mat]['color'],
                       markersize=10, capsize=4, capthick=2,
                       elinewidth=1.5, linestyle='none')
        # Configuración del subplot
        ax.set_title(f'{config_materiales[ref_mat]["label"]}', 
                    fontsize=18)#, pad=12, weight='bold')
        ax.grid(True, linestyle=':', alpha=0.6)
        #ax.axhline(y=3, color='red', linestyle='--', alpha=0.5, linewidth=1)
        #ax.set_ylim(-0.5, 10)
        ax.set_ylim(0, y_max)  # Mismo rango Y para todos
        ax.set_xlim(40, 260)
        ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Configuración de ejes
    for ax in axs[:, 0]:
        ax.set_ylabel(r'$\sigma$', fontsize=16, labelpad=10)
    for ax in axs[2, :]:
        ax.set_xlabel('Width (cm)', fontsize=16, labelpad=10)
    

    fig.legend(
    handles=[Line2D([0], [0], marker=config_materiales[mat]['marker'],
               color=config_materiales[mat]['color'],
               label=config_materiales[mat]['label'],
               markersize=12, linestyle='none')
             for mat in materiales],
    loc='center right',
    bbox_to_anchor=(1.1, 0.9),
    borderaxespad=0.,
    fontsize=18,
    )
    
    # Título general
    #plt.suptitle('Comparación de Significancias Estadísticas entre Materiales', 
                #fontsize=20, y=0.98, weight='bold')
    plt.tight_layout()
    plt.show()

# Ejecutar
graficar_comparaciones_significativas()
#import sys
#sys.exit()
# --------- Función para calcular tiempos y errores ----------
def obtener_conteos_sin_cambioDir(L, Mat):
    file_material = f'Flux-Mu-60min-{L}-{Mat}-Cor/'
    N_material = np.loadtxt(Carpet + file_material + 'Coinc_All.txt')
    return N_material

def obtener_conteos(L, Mat):
    file_material = f'Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    file_aire = f'Flux-Mu-60min-{L}-Air-UnDir-COR/'
    N_aire = np.loadtxt(Carpet + file_aire + 'Coinc_All.txt')
    N_material = np.loadtxt(Carpet + file_material + 'Coinc_All.txt')
    return N_aire, N_material

def calcular_tiempo_observacion(sigmas, sigma_errs, espesores, Mat_ref="Air"):
    t0_values = []
    for L in espesores:
        N_aire, _ = obtener_conteos(L, Mat_ref)
        N_sindir = len(obtener_conteos_sin_cambioDir(L, Mat_ref))
        t0 = len(N_aire) / N_sindir if N_sindir > 0 else 0
        t0_values.append(t0)
    t0_values = np.array(t0_values)
    tiempos = 9 * t0_values / (sigmas**2)
    tiempo_errs = 9 * (2 * t0_values * sigma_errs) / (sigmas**3)
    return tiempos, tiempo_errs, t0_values

def graficar_comparacion_tiempo_observacion(umbral_significancia=3.0):
    fig, axs = plt.subplots(3, 2, figsize=(20, 17), dpi=100)
    fig.subplots_adjust(hspace=0.3, wspace=0.2, right=0.85)

    config_materiales = {
        'Al': {'marker': 'o', 'color': 'black', 'label': 'Al'},
        'Con': {'marker': 's', 'color': 'orange', 'label': 'Concrete'},
        'Air': {'marker': 'D', 'color': 'red', 'label': 'Air'},
        'Wa': {'marker': '^', 'color': 'blue', 'label': 'Water'},
        'Fe': {'marker': 'v', 'color': 'green', 'label': 'Fe'},
        'Pb': {'marker': 'P', 'color': 'grey', 'label': 'Pb'}}

    factor_amplificacion_error = 1

    def calcular_tiempo_con_error(ref_sigma, ref_err, sigma, err):
        with np.errstate(divide='ignore', invalid='ignore'):
            diff = np.abs(ref_sigma - sigma)
            if diff == 0 or np.isnan(diff):
                return np.nan, np.nan
            sigma_comb = diff / np.sqrt(ref_err**2 + err**2)
            if sigma_comb <= 0:
                return np.nan, np.nan
            tiempo = (umbral_significancia / sigma_comb) ** 2
            dt_rel = 2 * (ref_err**2 + err**2) / (diff**2)
            dt = tiempo * np.sqrt(dt_rel)
            return tiempo, dt * factor_amplificacion_error

    # Determinar t_max para ajustar escala
    t_max = 0
    for ref_mat in materiales:
        ref_sigmas, ref_errs = results[ref_mat]
        for mat in materiales:
            if mat == ref_mat:
                continue
            sigmas, errs = results[mat]
            for i in range(len(espesores)):
                try:
                    t, _ = calcular_tiempo_con_error(ref_sigmas[i], ref_errs[i],
                                                     sigmas[i], errs[i])
                    if not np.isnan(t) and t > t_max:
                        t_max = t
                except:
                    continue
    t_max = min(1, t_max * 1.2)

    for idx, ref_mat in enumerate(materiales):
        ax = axs[idx // 2, idx % 2]
        ref_sigmas, ref_errs = results[ref_mat]

        for mat in materiales:
            if mat == ref_mat:
                continue
            sigmas, errs = results[mat]
            tiempos, tiempos_err = [], []
            for i in range(len(espesores)):
                try:
                    t, dt = calcular_tiempo_con_error(ref_sigmas[i], ref_errs[i],
                                                      sigmas[i], errs[i])
                    tiempos.append(t)
                    tiempos_err.append(dt)
                except:
                    tiempos.append(np.nan)
                    tiempos_err.append(np.nan)
            ax.errorbar(espesores, tiempos, yerr=tiempos_err,
                        fmt=config_materiales[mat]['marker'],
                        color=config_materiales[mat]['color'],
                        markersize=10, capsize=4, capthick=2,
                        elinewidth=1.5, linestyle='none')

        ax.set_title(f'{config_materiales[ref_mat]["label"]}', fontsize=18)
        ax.set_yscale('log')
        ax.set_xlim(40, 260)
        ax.set_ylim(0.001, 10)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.tick_params(axis='both', which='major', labelsize=12)

    for ax in axs[:, 0]:
        ax.set_ylabel('Observation time (a.u.)', fontsize=16)
    for ax in axs[2, :]:
        ax.set_xlabel('Width (cm)', fontsize=16)

    fig.legend(
        handles=[Line2D([0], [0], marker=config_materiales[mat]['marker'],
                        color=config_materiales[mat]['color'],
                        label=config_materiales[mat]['label'],
                        markersize=12, linestyle='none')
                 for mat in materiales],
        loc='center right',
        bbox_to_anchor=(1.1, 0.9),
        borderaxespad=0.,
        fontsize=18,
    )

    plt.tight_layout()
    plt.show()

# Ejecutar
graficar_comparacion_tiempo_observacion()

# Umbral de significancia
umbral_significancia = 3.0

# Espesor de interés
espesor_objetivo = 250
idx_250 = 4

# Obtener sigma y errores de cada material en ese espesor
sigma_Al, err_Al = results["Fe"][0][idx_250], results["Fe"][1][idx_250]
sigma_Pb, err_Pb = results["Pb"][0][idx_250], results["Pb"][1][idx_250]

# Cálculo del tiempo de observación
def calcular_tiempo_con_error(ref_sigma, ref_err, sigma, err):
    with np.errstate(divide='ignore', invalid='ignore'):
        diff = np.abs(ref_sigma - sigma)
        if diff == 0 or np.isnan(diff):
            return np.nan, np.nan
        sigma_comb = diff / np.sqrt(ref_err**2 + err**2)
        if sigma_comb <= 0:
            return np.nan, np.nan
        tiempo = (umbral_significancia / sigma_comb) ** 2
        dt_rel = 2 * (ref_err**2 + err**2) / (diff**2)
        dt = tiempo * np.sqrt(dt_rel)
        return tiempo, dt

# Ejecutar cálculo
tiempo_obs, error_obs = calcular_tiempo_con_error(sigma_Al, err_Al, sigma_Pb, err_Pb)

# Mostrar resultado
print(f"Tiempo de observación necesario para distinguir Al y Pb a {espesor_objetivo} cm:")
print(f"  Tiempo = {tiempo_obs:.2f} ± {error_obs:.5f} horas")'''