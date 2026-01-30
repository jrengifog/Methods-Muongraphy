import numpy as np
import matplotlib.pyplot as plt

# ===== Datos de materiales =====
materiales = [
    {"nombre": "Air",       "Z": 7.3,  "A": 14.6,   "rho": 0.001225},
    {"nombre": "Water",     "Z": 7.42, "A": 18.015, "rho": 1.0},
    {"nombre": "Al",        "Z": 13,   "A": 26.98,  "rho": 2.70},
    {"nombre": "Concrete",  "Z": 11,   "A": 22,     "rho": 2.3},
    {"nombre": "Fe",        "Z": 26,   "A": 55.85,  "rho": 7.87},
    {"nombre": "Pb",        "Z": 82,   "A": 207.2,  "rho": 11.34}
]

# ===== Parámetros =====
E_mu = 11000  # MeV (11 GeV)
p_mu = E_mu  # MeV/c, aproximación relativista
beta = 1.0
distancia_cm = 250  # longitud en cm

# ===== Función para X0 =====
def X0_radiacion(Z, A):
    return (716.4 * A) / (Z * (Z + 1) * np.log(287 / np.sqrt(Z)))  # g/cm²

# ===== Calcular theta para cada material =====
Z_vals, theta_vals = [], []

for m in materiales:
    Z = m["Z"]
    A = m["A"]
    rho = m["rho"]

    X0 = X0_radiacion(Z, A)
    X_gcm2 = distancia_cm * rho
    theta = (19.2 / (beta * p_mu)) * np.sqrt(X_gcm2 / X0) * (1 + 0.038 * np.log(X_gcm2 / X0))
    theta_deg = np.degrees(theta)  # convertir a grados

    Z_vals.append(Z)
    theta_vals.append(theta_deg)

# ===== Graficar =====
plt.figure(dpi=300)
plt.plot(Z_vals, theta_vals, 'o-', color='red', markerfacecolor='red',
         markeredgewidth=0.8, markersize=5, linewidth=1)
for i, m in enumerate(materiales):
    plt.text(Z_vals[i] + 0.5, theta_vals[i], m["nombre"], fontsize=8)

plt.xlabel("$Z$", fontsize=14)
plt.ylabel(r"$\theta$ (Degree)", fontsize=14)
plt.grid(True, which="both", ls="--", lw=0.5)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.show()