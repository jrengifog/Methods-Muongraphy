import numpy as np
import matplotlib.pyplot as plt

# === Fuente para graficar ===
font = {
    'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 16,
}

# === Ruta base ===
File = 'D:/fisica pc/documentos/Carpeta-share/dEdX_R/'

# === Leer archivos ===
DataAir = np.loadtxt(File + 'muE_air_dry_1_atm.txt')
DataPb = np.loadtxt(File + 'muE_lead_Pb.txt')
DataPol = np.loadtxt(File + 'muE_polyethylene.txt')
DataFe = np.loadtxt(File + 'muE_iron_Fe.txt')
DataAl = np.loadtxt(File + 'muE_aluminum_Al.txt')
DataCon = np.loadtxt(File + 'muE_shielding_concrete.txt')
DataWa = np.loadtxt(File + 'muE_water_liquid.txt')

# === Función para extraer energía y pérdida total ===
def D_mu_sh(Data):
    E = Data[:, 0]
    dEdX = Data[:, 7]
    return E, dEdX

# === Obtener datos ===
E_Air, dEdX_Air = D_mu_sh(DataAir)
E_Pb, dEdX_Pb = D_mu_sh(DataPb)
E_Pol, dEdX_Pol = D_mu_sh(DataPol)
E_Fe, dEdX_Fe = D_mu_sh(DataFe)
E_Al, dEdX_Al = D_mu_sh(DataAl)
E_Con, dEdX_Con = D_mu_sh(DataCon)
E_Wa, dEdX_Wa = D_mu_sh(DataWa)

# === Gráfico 1: Comparación entre materiales ===
plt.figure(figsize=(7, 5), dpi=300)
plt.plot(E_Wa, dEdX_Wa, marker='^', linestyle='None',  label="$H_2O$", color = 'blue',markersize=5)
plt.plot(E_Pol, dEdX_Pol, marker='D', linestyle='None', label="$[C_6H_5CHCH_2]_n$", color = 'brown',markersize=3)
plt.plot(E_Con, dEdX_Con, marker='s', linestyle='None', label="Concrete", color = 'orange',markersize=5)
plt.plot(E_Al, dEdX_Al, marker='o', linestyle='None', label="Al", color = 'black',markersize=5)
plt.plot(E_Fe, dEdX_Fe, marker='v', linestyle='None', label="Fe", color = 'green',markersize=5)
plt.plot(E_Pb, dEdX_Pb, marker='P', linestyle='None', label="Pb", color = 'gray',markersize=5)

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$E_\mu$ (MeV)", fontdict=font)
plt.ylabel(r"dE/dX (MeV cm$^2$ g$^{-1}$)", fontdict=font)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper left', fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()

# === Función para extraer componentes ===
def D_mu_sh_components(Data):
    E = Data[:, 0]
    dEdXT = Data[:, 7]
    dEdXI = Data[:, 2]
    dEdXB = Data[:, 3]
    dEdXPP = Data[:, 4]
    dEdXPN = Data[:, 5]
    dEdXRL = Data[:, 6]
    return E, dEdXT, dEdXI, dEdXB, dEdXPP, dEdXPN, dEdXRL

# === Obtener componentes de pérdida para polietileno ===
E_Pol, dEdX_T, dEdX_I, dEdX_B, dEdX_PP, dEdX_PN, dEdX_RL = D_mu_sh_components(DataPol)

# === Gráfico 2: Componentes en Polietileno ===
plt.figure(figsize=(7, 5), dpi=300)
plt.plot(E_Pol, dEdX_T, marker='o', linestyle='None', label="Total",markersize=5)
plt.plot(E_Pol, dEdX_I, marker='s', linestyle='None', label="Ionization",markersize=5)
plt.plot(E_Pol, dEdX_B, marker='^', linestyle='None', label="Bremsstrahlung",markersize=5)
plt.plot(E_Pol, dEdX_PP, marker='D', linestyle='None', label="Pair-Production",markersize=5)
plt.plot(E_Pol, dEdX_PN, marker='v', linestyle='None', label="Photo-Nuclear",markersize=5)
plt.plot(E_Pol, dEdX_RL, marker='>', linestyle='None', label="Radiative Loss",markersize=5)

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$E_\mu$ (MeV)", fontdict=font)
plt.ylabel(r"dE/dX (MeV cm$^2$ g$^{-1}$)", fontdict=font)
plt.ylim([1, 20000])
plt.xlim([1e3, 1e8])
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper left', fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()