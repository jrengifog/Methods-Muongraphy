import numpy as np
import matplotlib.pyplot as plt

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

File = 'D:/fisica pc/documentos/Carpeta-share/dEdX_R/'
DataAir = np.loadtxt(File + 'muE_air_dry_1_atm.txt')
DataPb = np.loadtxt(File + 'muE_lead_Pb.txt')
DataPol = np.loadtxt(File + 'muE_polyethylene.txt')
DataFe = np.loadtxt(File + 'muE_iron_Fe.txt')
DataAl = np.loadtxt(File + 'muE_aluminum_Al.txt')
DataCon = np.loadtxt(File + 'muE_shielding_concrete.txt')
DataWa = np.loadtxt(File + 'muE_water_liquid.txt')

def D_mu_sh(Data):
    E = []
    dEdX = []
    for i in range(0, len(Data)):
        E.append(Data[i, 0])
        dEdX.append(Data[i, 7])
    
    return E,dEdX

E_Air, dEdX_Air = D_mu_sh(DataAir)
E_Pb, dEdX_Pb = D_mu_sh(DataPb)
E_Pol, dEdX_Pol = D_mu_sh(DataPol)
E_Fe, dEdX_Fe = D_mu_sh(DataFe)
E_Al, dEdX_Al = D_mu_sh(DataAl)
E_Con, dEdX_Con = D_mu_sh(DataCon)
E_Wa, dEdX_Wa = D_mu_sh(DataWa)

plt.figure(dpi=300)
plt.plot(E_Wa, dEdX_Wa, 'blue', label="H2O")
plt.plot(E_Pol, dEdX_Pol,'green', label="$[C_6H_5CHCH_2]_n$")
plt.plot(E_Con, dEdX_Con, 'indigo', label="Concrete")
#plt.plot(E_Air, dEdX_Air,'blue', label="Air")
plt.plot(E_Al, dEdX_Al, 'red', label="Al")
plt.plot(E_Fe, dEdX_Fe,'chocolate', label="Fe")
plt.plot(E_Pb, dEdX_Pb,'grey', label="Pb")

#plt.ylim([10, 20000])
#plt.xlim([400000, 1000000000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel("$E_\mu$ (MeV)", fontdict=font)
plt.ylabel("dE/dX (MeV $cm^2$ $g^{-1}$)", fontdict=font)
plt.legend(loc='upper left', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

File = 'D:/fisica pc/documentos/Carpeta-share/dEdX_R/'
DataPol = np.loadtxt(File + 'muE_polyethylene.txt')
#DataPol_2 = np.loadtxt('C:/Users/Fisica/Documents/Carpeta-share/dEdX_Muon_polystyrene.txt')

def D_mu_sh(Data):
    E = []
    dEdXT,dEdXI,dEdXB, dEdXPP, dEdXPN, dEdXRL  = [],[],[],[],[],[]
    for i in range(0, len(Data)):
        E.append(Data[i, 0])
        dEdXT.append(Data[i, 7])
        dEdXI.append(Data[i, 2])
        dEdXB.append(Data[i, 3])
        dEdXPP.append(Data[i, 4])
        dEdXPN.append(Data[i, 5])
        dEdXRL.append(Data[i, 6])
        
    
    return E,dEdXT,dEdXI,dEdXB,dEdXPP,dEdXPN,dEdXRL 

E_Pol, dEdX_T, dEdX_I, dEdX_B, dEdX_PP, dEdX_PN, dEdX_RL  = D_mu_sh(DataPol)
#E_Pol2, dEdX_Pol_T2, dEdX_I2 = D_mu_sh(DataPol_2)


plt.figure()
#plt.plot(E_Wa, dEdX_Wa, 'blue', label="H2O")
plt.plot(E_Pol, dEdX_T,'green', label="Total")
plt.plot(E_Pol, dEdX_I,'b--', label="Ionization")
plt.plot(E_Pol, dEdX_B,'brown', label="Bremsstrahlung")
plt.plot(E_Pol, dEdX_PP,'purple', label="Pair-Production")
plt.plot(E_Pol, dEdX_PN,'red', label="Photo-Nuclear")
plt.plot(E_Pol, dEdX_RL,'orange', label="Radiation Loss")




plt.ylim([1, 20000])
plt.xlim([1000, 100000000])
#plt.plot(E_Con, dEdX_Con, 'indigo', label="Concrete")
#plt.plot(E_Air, dEdX_Air,'blue', label="Air")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("$E\mu$ (MeV)", fontdict=font)
plt.ylabel("dE/dX (MeV $cm^2$ $g^{-1}$)", fontdict=font)

plt.legend(loc='upper left', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
#plt.savefig('D:/Documentos Javier 2022/dEdX-distContr.png', dpi=300)
exit()

#
# Define constants
Z = 7.0  # atomic number of polystyrene
A = 14.0  # atomic mass of polystyrene
I = 82.0  # mean excitation energy of polystyrene, in eV
m_mu = 105.7  # muon mass, in MeV/c^2
density = 1.05  # density of polystyrene, in g/cm^3
re = 2.81794032e-13  # classical electron radius, in cm
alpha = 1/137  # fine structure constant
hbar = 6.582119514e-22  # Planck constant divided by 2*pi, in MeV*s

m_e =  9.10938356 * 10**(-28)# kg
c = 299792458 * 10**(2) #cm/s

# Define energy range
E = np.logspace(2, 5, 100)

# Calculate contributions from different scattering and energy loss mechanisms
bethe_bloch = np.zeros_like(E)
density_correction = np.zeros_like(E)
radiative_energy_loss = np.zeros_like(E)

for i in range(len(E)):
    gamma = (E[i] + m_mu) / m_mu
    beta = np.sqrt(1 - 1 / gamma ** 2)

    W_max = 2 * m_e * beta ** 2 * gamma ** 2 / (1 + 2 * gamma * m_mu / m_e + (m_mu / m_e)**2)
    ln_term = np.log(2 * m_e * beta ** 2 * gamma ** 2 * W_max * I ** 2 / (I ** 2 - (m_e * c ** 2) ** 2))
    delta = 0
    if beta >= 0.1:
        delta = 4.6052 + np.log(beta ** 2 * gamma ** 2) - np.log(I / (m_e * c ** 2))
    K = 4 * np.pi * re ** 2 * m_e * c ** 2 * Z / A
    S = density * K * Z / A * (1 / beta ** 2 - 1) * (0.5 * ln_term - beta ** 2 - 0.5 * delta)
    bethe_bloch[i] = S

# Calculate density correction term
    X = np.log10(beta * gamma)
    C = -5 + 0.5 * X - 0.5 * X ** 2 + 0.33 * X ** 3 - 0.14 * X ** 4
    density_correction[i] = C

    # Calculate radiative energy loss
    tau = 5.2e-17  # muon lifetime, in s
    beta_gamma = beta * gamma
    T_max = 2 * m_e * beta_gamma ** 2 / (1 + 2 * gamma * m_mu / m_e + (m_mu / m_e) ** 2)
    eta = T_max / (tau * beta_gamma ** 2)
    rad_energy_loss = -2 * m_e * c ** 2 * (beta ** 2 * gamma ** 2) * (np.log(eta) - beta ** 2)
    radiative_energy_loss[i] = rad_energy_loss

fig1 = plt.figure()
# Plot contributions from different scattering and energy loss mechanisms
plt.plot(E, bethe_bloch, label='Bethe-Bloch')
plt.plot(E, density_correction, label='Density Correction')
plt.plot(E, radiative_energy_loss, label='Radiative Energy Loss')
plt.plot(E, bethe_bloch + density_correction + radiative_energy_loss, label='Total')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Muon Energy (MeV)', fontdict=font)
plt.ylabel('Energy Loss (MeV/cm)', fontdict=font)

plt.legend(loc='upper left', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)



def energy_loss(muon_energy):
    # Constants for polystyrene
    ionization_constant = 2.3  # MeV/cm
    bremsstrahlung_constant = 0.024  # MeV/cm
    nuclear_interaction_constant = 0.003  # MeV/cm

    # Calculate energy loss contributions
    ionization_loss = ionization_constant * muon_energy
    bremsstrahlung_loss = bremsstrahlung_constant * np.log(muon_energy)
    nuclear_interaction_loss = nuclear_interaction_constant * muon_energy**0.75

    # Total energy loss
    total_loss = ionization_loss + bremsstrahlung_loss + nuclear_interaction_loss

    return ionization_loss, bremsstrahlung_loss, nuclear_interaction_loss, total_loss

# Generate muon energy values from 1 MeV to 1000 MeV
muon_energy_values = np.logspace(1,5, num=100)

# Calculate energy loss for each muon energy value
ionization_losses = []
bremsstrahlung_losses = []
nuclear_interaction_losses = []
total_losses = []

for energy in muon_energy_values:
    ionization, bremsstrahlung, nuclear_interaction, total = energy_loss(energy)
    ionization_losses.append(ionization)
    bremsstrahlung_losses.append(bremsstrahlung)
    nuclear_interaction_losses.append(nuclear_interaction)
    total_losses.append(total)

fig2 = plt.figure()
# Create the plot
plt.plot(muon_energy_values, ionization_losses, label='Ionization')
plt.plot(muon_energy_values, bremsstrahlung_losses, label='Bremsstrahlung')
plt.plot(muon_energy_values, nuclear_interaction_losses, label='Nuclear Interaction')
plt.plot(muon_energy_values, total_losses, label='Total Loss')

plt.xlabel('Muon Energy (MeV)', fontdict=font)
plt.ylabel('Energy Loss (MeV/cm)', fontdict=font)
plt.title('Muon Energy Loss in Polystyrene')

plt.xscale('log')
plt.yscale('log')

plt.legend(loc='upper left', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.grid(True)
