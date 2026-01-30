import numpy as np
import matplotlib.pyplot as plt

#Fuente para graficar
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

# Define the constants
area = 2500 * 1e-4 # in m^2
time = 1 * 3600 # in s
init_particles = (2 * np.pi * area * time)

# Define the primary cosmic ray nuclei spectra data
energy = np.logspace(2, 6, 1000) # in eV
proton_flux = 0.1 * ((energy/1e9)**-2.7) * init_particles # in particles
helium_flux = 0.0045 * ((energy/1e9)**-2.65) * init_particles # in particles
carbon_flux = 0.000066 * ((energy/1e9)**-2.75) * init_particles # in particles
oxygen_flux = 0.000040 * ((energy/1e9)**-2.7) * init_particles # in particles
iron_flux = 0.0000027 * ((energy/1e9)**-2.6) * init_particles # in particles
lead_flux = 0.0000003 * ((energy/1e9)**-2.4) * init_particles # in particles


# Configuración de figura
plt.figure(dpi=300, figsize=(7, 5))

# Lista de núcleos y sus estilos
nucleos = [
    {'label': 'p',  'flux': proton_flux, 'color': 'blue',    'marker': 'o'},
    {'label': 'He', 'flux': helium_flux, 'color': 'green',   'marker': 's'},
    {'label': 'C',  'flux': carbon_flux, 'color': 'orange',  'marker': '^'},
    {'label': 'O',  'flux': oxygen_flux, 'color': 'red',     'marker': 'v'},
    {'label': 'Fe', 'flux': iron_flux,   'color': 'purple',  'marker': 'd'},
    {'label': 'Pb', 'flux': lead_flux,   'color': 'brown',   'marker': 'P'},
]

# Gráfica log-log sin líneas, solo marcadores
#for nucleo in nucleos:
#    plt.loglog(
#        energy, nucleo['flux'],
#        marker=nucleo['marker'],
#        color=nucleo['color'],
#        linestyle='None',
#        markersize=4,
#        label=nucleo['label']
#    )
skip = 25  # muestra solo 1 de cada 25 puntos
for nucleo in nucleos:
    plt.loglog(
        energy[::skip], nucleo['flux'][::skip],
        marker=nucleo['marker'],
        color=nucleo['color'],
        linestyle='None',
        markersize=6,
        label=nucleo['label']
    )

plt.xlabel("Energy (GeV)", fontdict=font)
plt.ylabel("Flux ((m$^2$ sr s GeV)$^{-1}$)", fontdict=font)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Leyenda fuera del gráfico
plt.legend(loc='upper right', fontsize=12)#bbox_to_anchor=(1.05, 1), 
plt.tight_layout()
plt.show()
#plt.savefig('D:/Documentos Javier 2022/Flux_primary_CR_0.png', dpi=300)
#plt.figure(dpi=300)
#plt.show()
exit()

AlT = 16100  # [cm] ObsLev
BX = 23.0048
BZ = -0.498
Rig = 12.17  # GV in Lima
Area = 50*50  # [cm2]
Time = 1*3600  # [s]
N0 = (2*np.pi*Area*Time)
#print(N0)
Nnucleos = 26
# H, He, Li, Be, B, C, N , O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe
Z = list(range(1, 27))  # np.arange(27)
# print(Z)
# exit()
COR_ID_F = [14, 402, 703, 904, 1105, 1206, 1407, 1608, 1909, 2010, 2311, 2412, 2713, 2814, 3115, 3216, 3517, 3919,
            4018, 4020, 4521, 4822, 5123, 5224, 5525, 5626]
MASS_F = [0.938272, 3.73338244805557, 6.53478032991106, 8.39429044902688, 10.2536061929541, 11.1787903246739,
          13.045071978869, 14.898326507629, 17.6899146520668, 18.6173579550734, 21.4080199431823, 22.3362803688324,
          25.1263356296296, 26.0553153433303, 28.8449660324983, 29.7745989328225, 32.5639816988633, 36.2834316370329,
          37.2107457840596, 37.2142385732562, 41.8605295331555, 44.6483661801865, 47.4401999906342, 48.3681334024753,
          51.1598095147594, 52.0885229269484]

A0_F = [1.151, 7.19, 2.08, 4.74, 8.95, 1.06, 2.35, 1.57, 3.28, 4.6, 7.54, 8.01, 1.15, 7.96, 2.7, 2.29, 2.94, 8.36,
        5.36, 1.47, 3.04, 1.13, 6.31, 1.36, 1.35, 1.78]
E0_F = [5, 6, 7, 8, 8, 6, 7, 6, 8, 7, 8, 7, 7, 7, 8, 7, 8, 8, 8, 7, 8, 7, 8, 7, 7, 6]
Index_F = [2.77, 2.64, 2.54, 2.75, 2.95, 2.66, 2.72, 2.68, 2.69, 2.64, 2.66, 2.64, 2.66, 2.75, 2.69, 2.55, 2.68, 2.64,
           2.65, 2.70, 2.64, 2.61, 2.63, 2.67, 2.46, 2.60]

Nshow = 0
# Flux: j(E) = j0 * E^(-gamma)
Llimit = 100      #GeV
userllimit = Llimit
Ulimit = 1000000  #GeV

NSH = []
for i in range(0, len(Z)):
    J0 = A0_F[i] * 10 ** (E0_F[i]*(-1))
    G0_F = Index_F[i]*(-1)
    A_F = G0_F + 1
    #print(J0)
    if (Rig):
        P0 = Z[i] * Rig
        Llimit = np.sqrt(P0 ** 2 + MASS_F[i] ** 2)
        #print(Llimit)
        if (Llimit < userllimit):
            Llimit = userllimit
    if (Llimit < MASS_F[i]):
        Llimit = (COR_ID_F[i]/100 - Z[i]) * 0.9396 + Z[i] * 0.9383
        #print(Llimit)
    Nshow = int(N0 * (J0 / A_F) * ((Ulimit/1000.)**A_F - (Llimit/1000.)**A_F)) + 1
    print(Nshow, COR_ID_F[i], Llimit, Index_F[i])
    NSH.append(Nshow)
print(np.sum(NSH))
#print('--------------------')

ENERGY = np.linespace(np.logspace(10,1000000, 100))#list(range(100, 100000))#np.logspace(2,6, 100)#
NSHOW_H = []
NSHOW_He = []
NSHOW_Li = []
NSHOW_Be = []
NSHOW_B = []
for j in range(0,len(ENERGY)):
    NSHOW1 =  N0 * ((A0_F[j] * 10 ** (E0_F[j]*(-1))) / Index_F[j]*(-1)) * ((ENERGY[j]/1000.)**(Index_F[j]*(-1)))
    #(NSH[0])*ENERGY[j]**(Index_F[0]*(-1))#(10 ** (E0_F[1]*(-1)))*ENERGY[j]**(Index_F[1]*(-1))#N0*A0_F[1] * 10 ** (E0_F[1]*(-1))*(ENERGY[j]/1000.)**(Index_F[1]*(-1)+1)  #((N0 * (A0_F[1] * 10 ** (E0_F[1]*(-1))) * ((ENERGY[j]/1000.)**(Index_F[1]*(-1)+1))) + 1)
    #NSHOW2 =  (NSH[1])*ENERGY[j]**(Index_F[1]*(-1))#(10 ** (E0_F[8]*(-1)))* ENERGY[j]**(Index_F[8]*(-1))#N0*A0_F[8] * 10 ** (E0_F[8]*(-1))*(ENERGY[j]/1000.)**(Index_F[8]*(-1)+1)  #NSHOW2 = ((N0 * (A0_F[2] * 10 ** (E0_F[2] * (-1))) * ((ENERGY[j]/1000.)**(Index_F[2]*(-1)+1))) + 1)
    #NSHOW3 = ((N0 * (A0_F[3] * 10 ** (E0_F[3] * (-1))) * ((ENERGY[j]/1000.)**(Index_F[3]*(-1)+1))) + 1)
    #NSHOW4 = ((N0 * (A0_F[4] * 10 ** (E0_F[4] * (-1))) * ((ENERGY[j]/1000.)**(Index_F[4]*(-1)+1))) + 1)
    #NSHOW5 = ((N0 * (A0_F[25] * 10 ** (E0_F[25] * (-1))) * ((ENERGY[j]/1000.)**(Index_F[25]*(-1)+1))) + 1)
    NSHOW_H.append(NSHOW1)
    #NSHOW_He.append(NSHOW2)
    #NSHOW_Li.append(NSHOW3)
    #NSHOW_Be.append(NSHOW4)
    #NSHOW_B.append(NSHOW5)
# Nshow = int(N0 * (J0 / A_F) * ((Ulimit/1000.)**A_F - (Llimit/1000.)**A_F)) + 1'''
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }
fig = plt.figure()
plt.plot(ENERGY, NSHOW_H, label="p")
#plt.plot(ENERGY, NSHOW_He, label="He")
#plt.loglog(energy, carbon_flux, label="C")
#plt.loglog(energy, oxygen_flux, label="O")

plt.xlabel("Energy (GeV)", fontdict=font)
#plt.ylabel("Flux ((m$^2$ sr s GeV)$^{-1}$)", fontdict=font)
#plt.plot(ENERGY, NSHOW_H)
#plt.plot(ENERGY, NSHOW_He)
#plt.plot(ENERGY, NSHOW_Li)
#plt.plot(ENERGY, NSHOW_Be)
#plt.plot(ENERGY, NSHOW_B)
#plt.yscale('log')
plt.xscale('log')
plt.legend(loc='upper right', fontsize=10)
plt.show()