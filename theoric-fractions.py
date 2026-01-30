import numpy as np
import matplotlib.pyplot as plt
#Programa que de los tracks, me de los que pasan por d1 - d2, luego los que pasan,
#no pasan, deberían pasar y pasan, no deberian pasar y no pasan por el detector.

def count_muons_D(D):
    muxsho = []
    #fm = []
    for i in range(0, len(D)):
        muxsho.append(D[i][0])
    (umu, cmu) = np.unique(muxsho, return_counts=True)
    fmu = np.asarray((umu, cmu)).T
    #fm. append(fmu)
    return fmu

#Flux-Mu-60min-100-Pb-UnDir

#Fuente para graficar
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

Carpet= 'D:/fisica pc/documentos/Carpeta-share/'
#Mu_60min_Al_250cm_UnDirCor - 'Flux-Mu-60min-%d-%s-Cor/'
def FracMatvsAir(L, Mat):
    Carpet2 ='Flux-Mu-60min-%d-%s-UnDir-COR/' % (L, Mat)#
    Carpet1 = 'Flux-Mu-60min-%d-Air-UnDir-COR/'% (L)
    ##
    FS = np.loadtxt(Carpet + Carpet1 + 'Coinc_All.txt')
    #FS = np.loadtxt(Carpet + Carpet2 + 'coincidences_Bef.txt')
    F4 = np.loadtxt(Carpet + Carpet2 + 'Coinc_All.txt')
    Frac = len(F4)/len(FS)
    Er_Frac = ((((np.sqrt(len(FS))) / len(FS)) + (((np.sqrt(len(F4))) / len(F4)))) * (Frac))
    return Frac, Er_Frac

def FracMatvsD4(L, Mat):
    Carpet2 = 'Flux-Mu-60min-%d-%s-UnDir-COR/'% (L, Mat)#
    #Carpet1 = 'Flux-Mu-60min-%d-Air-UnDir-COR/'% (L)
    #'Flux-Mu-60min-%d-%s-UnDir-COR/'
    #FS = np.loadtxt(Carpet + Carpet1 + 'Coinc_All.txt')
    FS = np.loadtxt(Carpet + Carpet2 + 'coincidences_Bef.txt')
    F4 = np.loadtxt(Carpet + Carpet2 + 'Coinc_All.txt')
    Frac = len(F4)/len(FS)
    Er_Frac = ((((np.sqrt(len(FS))) / len(FS)) + (((np.sqrt(len(F4))) / len(F4)))) * (Frac))
    return Frac, Er_Frac

M_Al = 'Al'
M_Con = 'Con'
M_W = 'Wa'
M_Air = 'Air'
M_Pb = 'Pb'
M_Fe = 'Fe'

W_Mat = [50, 100,150,200,250]#[150,200,250]#

def FraqD4(Mat):
    FracM, ErF = [],[]
    for i in [50,100,150,200,250]:
        Frac, er = FracMatvsD4(i, Mat)
        FracM.append(Frac)
        ErF.append(er)
    return FracM, ErF

def FraqAir(Mat):
    FracM, ErF = [],[]
    for i in [50,100,150,200,250]:
        Frac, er = FracMatvsAir(i, Mat)
        FracM.append(Frac)
        ErF.append(er)
    return FracM, ErF

FracP, erP = FraqD4(M_Pb)



#SPb, Er = arr_Sig(M_Pb)


W_Mat = [50, 100, 150, 200, 250]

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111)
plt.plot(W_Mat, FracP, color='grey', label="Simulated Pb")
plt.errorbar(W_Mat, FracP, yerr =erP,fmt='.',ecolor = 'grey',color='grey')
#,marker='o', linestyle='dotted'
#def calculate_moliere_radius(energy, Z, A, D):
    # Fórmula generalizada para calcular el radio de Molière
#    Xo = 716.4* A/(D*(Z*(Z+1))*((ln(287/np.sqrt(Z)))))
#    Rm = (13.6 * 1e3/energy)*np.sqrt(1/Xo)*(1 + 0.038*ln(1/Xo))#(15/energy)*np.sqrt(1/Xo)#

   
    #Rm = (21 * 1e6) / (energy * atomic_number * (1 + 0.038 * np.log(energy * atomic_number)))
#    return Rm#

#def Ang_Theoria(energy, Z, A, D, X):
    # Fórmula generalizada para calcular el radio de Molière
#    Xo = 716.4* A/(D*(Z*(Z+1))*((ln(287/np.sqrt(Z)))))
#    Th = (13.6 * 10**(-3)/energy)*np.sqrt(1/(1*Xo))*(1 + 0.038*ln(1/(1*Xo)))#(15/energy)*np.sqrt(1/Xo)#

   
    #Rm = (21 * 1e6) / (energy * atomic_number * (1 + 0.038 * np.log(energy * atomic_number)))
#    return np.degrees(Th)

grosor_plomo = np.linspace(50, 250, 50)

# Supongamos una tasa inicial constante de muones
tasa_muones_inicial = 0.62#0.55#0.62  # Por ejemplo, 1000 muones por segundo

# Supongamos un coeficiente de absorción ficticio (solo para demostración)
# Esta es una fórmula simplificada y ficticia. Necesitarías usar una fórmula real basada en datos experimentales.
coef_absorcion = 0.0016#0.001#0.0013 #0.6 * np.exp(0.001 * grosor_plomo)

# Calcula la tasa de muones después de la absorción en función del grosor del plomo
tasa_muones_despues_absorcion = tasa_muones_inicial * np.exp(-coef_absorcion * grosor_plomo)

# Crear la gráfica
plt.plot(grosor_plomo, tasa_muones_despues_absorcion, label='Model', color='crimson', linestyle='dashed', linewidth=2)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.legend(loc='upper right')
#fig.legend(fontsize=10)
#fig.tight_layout()
#fig.subplots_adjust(right=0.80)
fig.text(0.5, 0.01, 'Width (cm)', ha='center', va='center', fontdict=font)
fig.text(0.01, 0.5, 'Nf/Ni', ha='center', va='center', rotation='vertical', fontdict=font)
