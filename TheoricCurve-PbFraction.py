import numpy as np
import matplotlib.pyplot as plt

# ---------------- Config ----------------
Carpet = 'D:/fisica pc/documentos/Carpeta-share/'

# funciones de IO (las mismas que usabas)
def FracMatvsD4(L, Mat):
    carpeta = f'Flux-Mu-60min-{L}-{Mat}-UnDir-COR/'
    FS = np.loadtxt(Carpet + carpeta + 'coincidences_Bef.txt')
    F4 = np.loadtxt(Carpet + carpeta + 'Coinc_All.txt')
    Frac = len(F4) / len(FS)
    Er_Frac = ((((np.sqrt(len(FS))) / len(FS)) +
                ((np.sqrt(len(F4))) / len(F4)))) * Frac
    return Frac, Er_Frac

def FraqD4(Mat, espesores):
    FracM, ErF = [], []
    for L in espesores:
        Frac, er = FracMatvsD4(L, Mat)
        FracM.append(Frac)
        ErF.append(er)
    return np.array(FracM), np.array(ErF)

# ---------------- Física: modelos ----------------
def exponential_model(L_cm, E0_MeV, dEdx_MeVcm2g, rho_gcm3):
    """Modelo exponencial simple usando lambda = E0/(dEdx * rho)."""
    lambda_cm = E0_MeV / (dEdx_MeVcm2g * rho_gcm3)
    return np.exp(-L_cm / lambda_cm)

def spectrum_cut_model(L_cm, dEdx_MeVcm2g, rho_gcm3,
                       E_lo_MeV=100., E_hi_MeV=1e7, n=2.6):
    """
    Fracción que sobrevive si el espectro es phi(E) ~ E^{-n}.
    E_min = dEdx * rho * L  (pérdida esperada en MeV)
    Calcula integral analítica de E^{-n}.
    """
    E_min = dEdx_MeVcm2g * rho_gcm3 * L_cm  # MeV
    E_min = max(E_min, E_lo_MeV*1e-12)  # evitar 0
    # integrales:
    if n == 1.0:
        # integral ~ ln(E)
        num = np.log(E_hi_MeV / E_min)
        den = np.log(E_hi_MeV / E_lo_MeV)
    else:
        num = (E_min**(1.0 - n) - E_hi_MeV**(1.0 - n)) / (1.0 - n)
        den = (E_lo_MeV**(1.0 - n) - E_hi_MeV**(1.0 - n)) / (1.0 - n)
    # si denominador 0 accidental, devolver 0
    if den == 0:
        return 0.0
    return num / den

# ---------------- Parámetros y datos ----------------
rho_Pb = 11.34          # g/cm3
E_mu_typ = 11000#6000.0       # MeV (3 GeV), para modelo exponencial
dEdx_typ = 3.126*11.35#3.0          # MeV cm^2 / g, aproximación MIP

W_Mat = np.array([50, 100, 150, 200, 250])  # cm
FracP, erP = FraqD4("Pb", W_Mat)

# modelo exponencial (como referencia)
Lgrid = np.linspace(50, 250, 300)
expo_curve = exponential_model(Lgrid, E_mu_typ, dEdx_typ, rho_Pb)

# ---------------- Ajuste del modelo espectral (n) por búsqueda sencilla ----------------
n_values = np.linspace(0.8, 1.3, 100)  # buscar n entre 1.5 y 4.0
best_n = None
best_sse = np.inf
best_curve = None

# límites de energía usados para la integral (MeV)
E_lo = 100.0      # energía mínima considerada en la población
E_hi = 1e6        # corte en alta energía (10 TeV -> 1e7 MeV)

for n in n_values:
    pred = np.array([spectrum_cut_model(L, dEdx_typ, rho_Pb, E_lo, E_hi, n) for L in W_Mat])
    sse = np.sum((pred - FracP)**2)
    if sse < best_sse:
        best_sse = sse
        best_n = n
        best_curve = np.array([spectrum_cut_model(L, dEdx_typ, rho_Pb, E_lo, E_hi, best_n) for L in Lgrid])

# imprime resultado ajuste
print(f"Best-fit n = {best_n:.3f}, SSE = {best_sse:.3e}")

# ---------------- Plot: datos, exponencial, modelo espectral ajustado ----------------
plt.figure(dpi=300, figsize=(7,5))

# datos Pb
plt.errorbar(W_Mat, FracP, yerr=erP, fmt='P', color='grey',
             markersize=6, linewidth=1.5, capsize=4, capthick=1.2,
             elinewidth=1.2, label='Pb (data)')

# curva exponencial simple
#plt.plot(Lgrid, expo_curve, 'r--', linewidth=2, label='Exponential Model (E=3GeV)')

# curva espectral (mejor n)
plt.plot(Lgrid, best_curve, color='tab:red', linestyle='-.', linewidth=2,
         label='Exponential model')

plt.xlabel('Width (cm)', fontsize=16)
plt.ylabel('Nf/Ni', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()
