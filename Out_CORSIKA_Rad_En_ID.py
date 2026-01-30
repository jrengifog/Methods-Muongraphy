import numpy as np
import matplotlib.pyplot as plt

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

File = 'D:/fisica pc/documentos/Carpeta-share/'
DataE2 = File + 'Proton_E2V0.shw.bz2'
DataE3 = File + 'Proton_E3V0.shw.bz2'
DataE4 = File + 'Proton_E4V0.shw.bz2'
DataE5 = File + 'Proton_E2-5V0.shw.bz2'
DataE6 = File + 'Proton_E2-6V0.shw.bz2' #'Flux60min_DMD.shw' #
#DataE7 = File + 'Proton_E2-6.shw.bz2'

def Read_Par(Data):
    X,Y,Z = [],[],[]
    CORSIKAID, shower_id = [],[]
    Rad_mu, E_mu, sh_mu, R= [], [], [], []
    Nmu = 0
    Masa_Mu = 0.106 #GeV/c2
    DRead = np.loadtxt(Data)
    
    for i in range(0, len(DRead)):
        #X.append(DRead[i, 4])
        #Y.append(DRead[i, 5])
        Z.append(DRead[i, 6])
        CORSIKAID.append(DRead[i,0])
        shower_id.append(DRead[i,7])
        R.append(np.sqrt(DRead[i, 4]**2+DRead[i, 5]**2))#+DRead[i, 6]**2)
        if DRead[i,0] == 5 or DRead[i,0] == 6:#Muons 
            X.append(DRead[i, 4]*0.00001)
            Y.append(DRead[i, 5]*0.00001)
        #PX.append(DRead[i, 1])
        #PY.append(DRead[i, 2])
        #PZ.append(DRead[i, 3])   
            #RAD = np.sqrt(DRead[i, 4]**2+DRead[i, 5]**2)
            RAD3D = np.sqrt((DRead[i, 4]*0.00001)**2+(DRead[i, 5]*0.00001)**2)#*0.00001#+DRead[i, 6]**2)
            MOMT = np.sqrt(DRead[i, 1]**2+DRead[i, 2]**2+DRead[i, 3]**2)
            Nmu = Nmu + 1 # Number of muons of all showers
            Rad_mu.append(RAD3D) # Muons Radii
            ENERGY_mu = np.sqrt(MOMT** 2 + Masa_Mu ** 2)
            E_mu.append(ENERGY_mu)
            sh_mu.append(DRead[i,7])
        
    return  CORSIKAID,shower_id,X,Y,Z, Nmu,Rad_mu, E_mu, sh_mu, R

#CorsikaId px py pz x y z shower_id prm_id prm_energy prm_theta prm_phi
C_ID_2, S_ID_2, X2,Y2,Z2, Nmu2,R_mu2, E_mu2, SH_mu2, RT2 = Read_Par(DataE2)
C_ID_3, S_ID_3, X3,Y3,Z3, Nmu3,R_mu3, E_mu3, SH_mu3, RT3 = Read_Par(DataE3)
C_ID_4, S_ID_4, X4,Y4,Z4, Nmu4,R_mu4, E_mu4, SH_mu4, RT4 = Read_Par(DataE4)
C_ID_5, S_ID_5, X5,Y5,Z5, Nmu5,R_mu5, E_mu5, SH_mu5, RT5 = Read_Par(DataE5)
C_ID_6, S_ID_6, X6,Y6,Z6, Nmu6,R_mu6, E_mu6, SH_mu6, RT6 = Read_Par(DataE6)

'''
plt.figure()
plt.plot(X6,Y6,'b.')
#plt.plot(X_E3,Y_E3,'r.')
#plt.plot(X_E4,Y_E4,'g.')
plt.xlabel("Positions(cm)")
plt.ylabel("Positions(cm)")
plt.xlim([-25000,25000])
plt.ylim([-25000,25000])
#print(ID_E2)
#print(ID_E3)
#print(ID_E4)

#plt.figure()
fig, ax = plt.subplots()
# syntax for 3-D projection
#ax = plt.axes()#projection ='3d')
# defining axes
z = R_mu6 #S_ID_6#R_mu6#Nmu6#Z6# 
x = X6
y = Y6
#Z = z.reshape(21, 21)
#ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
#ax.scatter_density(x, y)
plt.scatter(x, y, c=z, cmap='gnuplot_r')#'CMRmap_r')
clb=plt.colorbar()
clb.ax.tick_params(labelsize=10) 
clb.ax.set_title('Radius(Km)',fontsize=8)
plt.xlabel("X (Km)", fontdict=font)
plt.ylabel("Y (Km)", fontdict=font)
#plt.savefig('D:/Documentos Javier 2022/PositionAllParticles_Lima60sec_102-106.png', dpi=300)
#plt.show()
#exit()
X1 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

squad = [ ' ' ,'$\gamma$','$e^+$','$e^-$',' ','$\mu^+$', '$\mu^-$', '$\pi$0','$\pi$+', '$\pi$-', '$K^0$', '$K^+$',  '$K^-$', '$n^0$', '$p^+$', '$p^-$']
kwargs1 = dict(histtype='step', alpha=0.75, density= False, bins=range(1,17,1), linewidth=1.5)
fig, ax1 = plt.subplots(1,1)
#ax1.hist(C_ID_2, color='b', label="0.1 TeV", **kwargs1)
#ax1.hist(C_ID_3, color='r', label="1 TeV", **kwargs1)
#ax1.hist(C_ID_4, color='g', label="10 TeV", **kwargs1)
#ax1.hist(C_ID_5, color='purple', label="0.1 - 100TeV", **kwargs1)
ax1.hist(C_ID_6, color='r', label="0.1 - 1000TeV", **kwargs1)
#plt.xticks(C_ID_6,squad)
ax1.set_xlabel("Corsika ID", fontdict=font)
ax1.set_ylabel("Counts", fontdict=font)
ax1.set_yscale("log")
#ax1.legend(loc='upper right', fontsize=10)
ax1.set_xticks(X1,fontsize=10)
#ax1.set_yticks(fontsize=10)
ax1.set_xticklabels(squad, minor=False, rotation=0)
#lab = [item.get_text() for item in ax1.get_xticklabels(X1)]
#lab = ['Photon','Positron','Electron',' ','Antimuon', 'Muon', '$\pi$0','$\pi$+', '$\pi$-', '$K$0', '$K$+',  '$K$-', 'Neutrón', 'Proton', 'Antiproton']
#ax1.set_xticklabels(lab, minor=False, rotation=45)
kwargs1 = dict(histtype='step', alpha=0.75, density= False, bins=range(1,17,1), linewidth=1.5)
'''
#plt.title("Primary Cosmic Ray Nuclei Spectra (Area={} cm^2, Observation Time={} s)".format(area*1e4, time))
#plt.xlabel("Energy (GeV)", fontdict=font)
#plt.ylabel("Flux ((m$^2$ sr s GeV)$^{-1}$)", fontdict=font)
#plt.legend(loc='upper right', fontsize=10)
#plt.xticks(fontsize=10)
#plt.yticks(fontsize=10)
#plt.savefig('D:/Documentos Javier 2022/COR-ID-Range-from-flux.png', dpi=300)

#handles, labels = ax.get_legend_handles_labels()
#new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]

#plt.legend(handles=new_handles, labels=labels)

#exit()
def selection_sort(x):
    for i in range(len(x)):
        swap = i + np.argmin(x[i:])
        (x[i], x[swap]) = (x[swap], x[i])
    Rs = 0.68*len(x)
    RSig = []
    for i in range(0, round(Rs)):
        RSig.append(x[i])
    return RSig

RS_E2 = selection_sort(R_mu2)
RS_E3 = selection_sort(R_mu3)
RS_E4 = selection_sort(R_mu4)
RS_E5 = selection_sort(R_mu5)
RS_E6 = selection_sort(R_mu6)



from matplotlib.lines import Line2D

kwargs = dict(histtype='step', alpha=0.75, density= False,  linewidth=1.5, bins = np.linspace(0, 4, num=60))#,bins=range(0,3,1))
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(R_mu2, color='b', label="0.1 TeV", **kwargs)
ax.hist(R_mu3, color='r', label="1 TeV", **kwargs)
ax.hist(R_mu4, color='g', label="10 TeV", **kwargs)
ax.hist(R_mu5, color='purple', label="0.1-100 TeV", **kwargs)
ax.hist(R_mu6, color='orange', label="0.1-1000 TeV", **kwargs)

plt.yscale('log')
plt.xlabel("Radii(km)", fontdict=font)
plt.ylabel("Counts", fontdict=font)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
handles, labels = ax.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
plt.legend(handles=new_handles, labels=labels, fontsize=10)
#plt.savefig('D:/Documentos Javier 2022/Rad_EnergyFix_AllPart.png', dpi=300)

fig = plt.figure()
ax1 = fig.add_subplot(111)
kwargsE = dict(histtype='step', alpha=0.75, density= False,  linewidth=1.5, bins=np.logspace(np.log10(0.1),np.log10(20000),100))#bins=range(0,100,1))

ax1.hist(E_mu2, color='b', label="0.1 TeV", **kwargsE)
ax1.hist(E_mu3, color='r', label="1 TeV", **kwargsE)
ax1.hist(E_mu4, color='g', label="10 TeV", **kwargsE)
ax1.hist(E_mu5, color='purple', label="0.1-100 TeV", **kwargsE)
ax1.hist(E_mu6, color='orange', label="0.1-1000 TeV", **kwargsE)
#plt.hist(ID_E3, density=False, alpha=0.75)
#plt.hist(ID_E4, density=False, alpha=0.75)
#plt.xlim([0, 10000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel("E$\mu$ (GeV)", fontdict=font)
plt.ylabel("Counts", fontdict=font)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
handles, labels = ax1.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
plt.legend(handles=new_handles, labels=labels, fontsize=10)
#plt.savefig('D:/Documentos Javier 2022/Emu_EnergyFix_AllPart.png', dpi=300)

kwargs = dict(histtype='step', alpha=0.75, density= False, bins=range(0,2300,50))
plt.figure()
plt.hist(RS_E2, color='b', label="0.1 TeV", **kwargs)
plt.hist(RS_E3, color='r', label="1 TeV", **kwargs)
plt.hist(RS_E4, color='g', label="10 TeV", **kwargs)
plt.hist(RS_E5, color='purple', label="0.1-100 TeV", **kwargs)
plt.hist(RS_E6, color='orange', label="0.1-1000 TeV", **kwargs)
#plt.hist(ID_E3, density=False, alpha=0.75)
#plt.hist(ID_E4, density=False, alpha=0.75)
#plt.xlim([0, 10000])
plt.yscale('log')
plt.xlabel("Radii(cm)", fontdict=font)
plt.ylabel("Counts", fontdict=font)
plt.legend(loc='upper right', fontsize=10)

#exit()

print('Nmuones')
print(Nmu2) 
print(Nmu3)
print(Nmu4)
print(Nmu5)
print(Nmu6)  
print('Mean Rad')
print(np.mean(R_mu2))
print(np.max(R_mu2)) 
print(np.mean(R_mu3))
print(np.max(R_mu3))
print(np.mean(R_mu4))
print(np.max(R_mu4))
print(np.mean(R_mu5))
print(np.mean(R_mu6))
print(np.max(R_mu6)) 
#exit()
print('Rad 1 sig')
print(len(RS_E2))
print(len(RS_E3))
print(len(RS_E4))
print(len(RS_E5))
print(len(RS_E6)) 

kwargs = dict(histtype='step', alpha=0.75, density= False, bins=np.logspace(np.log10(0.1),np.log10(1000),100))#bins=range(0,100,1))
plt.figure()
plt.hist(E_mu2, color='b', label="0.1 TeV", **kwargs)
plt.hist(E_mu3, color='r', label="1 TeV", **kwargs)
plt.hist(E_mu4, color='g', label="10 TeV", **kwargs)
plt.hist(E_mu5, color='purple', label="0.1-100 TeV", **kwargs)
plt.hist(E_mu6, color='orange', label="0.1-1000 TeV", **kwargs)
#plt.hist(ID_E3, density=False, alpha=0.75)
#plt.hist(ID_E4, density=False, alpha=0.75)
#plt.xlim([0, 10000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel("E$\mu$ (GeV)", fontdict=font)
plt.ylabel("Counts", fontdict=font)
plt.legend(loc='upper right', fontsize=10)
print('E Mean')
print(np.mean(E_mu2)) 
print(np.mean(E_mu3))
print(np.mean(E_mu4))
print(np.mean(E_mu6))   
#exit()        

def Nmu_Rad_En_x_Show(SH_mu,Rad_mu,E_mu):
    (u,c) = np.unique(SH_mu, return_counts=True)
    fmu = np.asarray((u,c)).T
    muxsho, R_mean_mu, E_mean_mu = [], [], []
    for i in range(0, len(fmu)):
        muxsho.append(fmu[i][1])#Cuantas particulas hay por cada lluvia
    for d in set(SH_mu):#Mean of Radio per each shower_ID
        R_mean_mu.append(np.mean([Rad_mu[i] for i in range(len(SH_mu)) if SH_mu[i] == d]))
        E_mean_mu.append(np.mean([E_mu[i] for i in range(len(SH_mu)) if SH_mu[i] == d]))
    return muxsho, R_mean_mu, E_mean_mu
    
MuxSho_E2, RxSho_E2, ExSho_E2 = Nmu_Rad_En_x_Show(SH_mu2,R_mu2,E_mu2)
MuxSho_E3, RxSho_E3, ExSho_E3 = Nmu_Rad_En_x_Show(SH_mu3,R_mu3,E_mu3)
MuxSho_E4, RxSho_E4, ExSho_E4 = Nmu_Rad_En_x_Show(SH_mu4,R_mu4,E_mu4)
MuxSho_E5, RxSho_E5, ExSho_E5 = Nmu_Rad_En_x_Show(SH_mu5,R_mu5,E_mu5)
MuxSho_E6, RxSho_E6, ExSho_E6 = Nmu_Rad_En_x_Show(SH_mu6,R_mu6,E_mu6)

#plt.figure()
#plt.hist(MuxSho_E2, bins=10, color = 'green', histtype='step', alpha=0.75, density= False)
#plt.xlabel("Number of muons per shower", fontdict=font)
#plt.ylabel("Counts", fontdict=font)
#plt.xlim([0, 100])
#plt.yscale("log")

print('Mu Mean x Sho')
print(np.mean(MuxSho_E2)) 
print(np.mean(MuxSho_E3))
print(np.mean(MuxSho_E4))
print(np.mean(MuxSho_E5))
print(np.mean(MuxSho_E6))

kwargs = dict(histtype='step', alpha=0.75, density= False, bins=range(0,30000,800))
plt.figure()
plt.hist(RxSho_E2, color='b', label="0.1 TeV", **kwargs)
plt.hist(RxSho_E3, color='r', label="1 TeV", **kwargs)
plt.hist(RxSho_E4, color='g', label="10 TeV", **kwargs)
#plt.hist(RxSho_E5, color='purple', label="0.1-100 TeV", **kwargs)
plt.hist(RxSho_E6, color='orange', label="0.1-1000 TeV", **kwargs)

#plt.xlim([0, 10000])
plt.yscale('log')
plt.xlabel("Radii(cm)", fontdict=font)
plt.ylabel("Counts", fontdict=font)
plt.legend(loc='upper right', fontsize=10)

print('R Mu Mean x Sho')
print(np.mean(RxSho_E2)) 
print(np.max(RxSho_E2)) 
print(np.mean(RxSho_E3))
print(np.max(RxSho_E3))
print(np.mean(RxSho_E4))
print(np.max(RxSho_E4))
print(np.mean(RxSho_E5))
print(np.mean(RxSho_E6))
print(np.max(RxSho_E6))

print('R Mu Mean x Sho 1 Sig')
print(np.mean(selection_sort(RxSho_E2))) 
print(np.max(selection_sort(RxSho_E2))) 
print(np.mean(selection_sort(RxSho_E3)))
print(np.max(selection_sort(RxSho_E3)))
print(np.mean(selection_sort(RxSho_E4)))
print(np.max(selection_sort(RxSho_E4)))
print(np.mean(selection_sort(RxSho_E5)))
print(np.mean(selection_sort(RxSho_E6)))
print(np.max(selection_sort(RxSho_E6)))

kwargs = dict(histtype='step', alpha=0.75, density= False, bins=np.logspace(np.log10(0.1),np.log10(1000),100))
plt.figure()
plt.hist(ExSho_E2, color='b', label="0.1 TeV", **kwargs)
plt.hist(ExSho_E3, color='r', label="1 TeV", **kwargs)
plt.hist(ExSho_E4, color='g', label="10 TeV", **kwargs)
#plt.hist(ExSho_E5, color='purple', label="0.1-100 TeV", **kwargs)
plt.hist(ExSho_E6, color='orange', label="0.1-1000 TeV", **kwargs)

#plt.xlim([0, 10000])
plt.xscale('log')
plt.yscale('log')
plt.xlabel("$E_\mu$ (GeV)", fontdict=font)
plt.ylabel("Counts", fontdict=font)
plt.legend(loc='upper left', fontsize=10)

print('E Mu Mean x Sho')
print(np.mean(ExSho_E2)) 
print(np.mean(ExSho_E3))
print(np.mean(ExSho_E4))
print(np.mean(ExSho_E5))
print(np.mean(ExSho_E6))

'''
fig2.tight_layout()
#exit()
#print(muxsho)
#print(np.mean(muxsho))
figmu = plt.figure()
plt.hist(muxsho, bins=100, color = 'green', histtype='step', alpha=0.75, density= False)
plt.xlabel("Number of muons per shower", fontdict=font)
plt.ylabel("Counts", fontdict=font)
#plt.xlim([0, 100])
plt.yscale("log")
'''
#plt.show()
#exit()