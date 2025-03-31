#%% VSM old ways - 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import os
from sklearn.metrics import r2_score 
from mlognormfit import fit3
from mvshtools import mvshtools as mt
import re
def lineal(x,m,n):
    return m*x+n
        
#%% Leo Archivos 
# patron Superparamagnetico 1
SP2 = np.loadtxt(os.path.join('VSM','Patron_SPM2a.txt'), skiprows=12)
H_SP2 = SP2[:, 0]  # Gauss
m_SP2 = SP2[:, 1]  # emu

factor_Flavio= 6.92/6.902 # debido a calibracion posterior
m_SP2*=factor_Flavio
# mass_capsula = 1.19897  # g
# mass_capsula_c_SP1 = 1.47020  # g
mass_muestra_SP2 = 0.17659 # g
C_SP2 = 18  #concentracion estimada en g/L = kg/m³
# Normalizo momento por masa de NP
m_SP2 /= mass_muestra_SP2  # emu/g


#%% Generar señales anhisteréticas
H_anhist_SP2, m_anhist_SP2 = mt.anhysteretic(H_SP2, m_SP2)

# Graficar señales anhisteréticas
fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
ax.plot(H_SP2, m_SP2, '.-', label='SP2')
ax.plot(H_anhist_SP2, m_anhist_SP2, '.-', label='SP2 anhisteretica')

for a in [ax]:
    a.legend(ncol=2)
    a.grid()
    a.set_ylabel('m (emu/g)')
plt.xlabel('H (G)')
plt.show()

#%% Realizo fits en ciclos CON contribucion diamag
fit_SP2 = fit3.session(H_anhist_SP2, m_anhist_SP2, fname='SP2', divbymass=False)
fit_SP2.fix('sig0')
fit_SP2.fix('mu0')
fit_SP2.free('dc')
fit_SP2.fit()
fit_SP2.update()
fit_SP2.free('sig0')
fit_SP2.free('mu0')
fit_SP2.set_yE_as('sep')
fit_SP2.fit()
fit_SP2.update()
fit_SP2.save()
fit_SP2.print_pars()
H_SP2_fit = fit_SP2.X
m_SP2_fit = fit_SP2.Y
m_SP2_sin_diamag = m_anhist_SP2 - lineal(H_anhist_SP2, fit_SP2.params['C'].value, fit_SP2.params['dc'].value)

#%% Graficar resultados eliminando comportamiento diamagnético
fig, ax = plt.subplots(nrows=1, figsize=(8, 6),constrained_layout=True)

ax.plot(H_SP2, m_SP2, '.-', label='SP2')
ax.plot(H_anhist_SP2, m_anhist_SP2, '.-', label='SP2 anhisteretica')
ax.plot(H_anhist_SP2, m_SP2_sin_diamag, '-', label='SP2 s/ diamag')
ax.plot(H_SP2_fit, m_SP2_fit, '-', label='SP2 fit')

for a in [ax]:
    a.legend(ncol=1)
    a.grid()
    a.set_ylabel('m (emu/g)')
plt.xlabel('H (G)')
plt.show()

#%% Salvo ciclo VSM

data = np.column_stack((H_SP2_fit, m_SP2_fit))

# Guarda la matriz en un archivo de texto
np.savetxt('archivo.txt', data, fmt=('%e','%e'),
           header='',comment=,delimiter='\t')
