'''
Write better code
'''

import cantera as ct
import numpy as np

spec = 'C6H3'
q = spec+':1'

P1 = ct.one_atm  # convert torr to Pa

T_range = [300, 6000]
stepsize = 100
n = (T_range[1]-T_range[0])/(stepsize)+1
T = np.append(298.15, np.linspace(T_range[0], T_range[1],n))

gas = ct.Solution('generated_mech.cti')

gas.TPX = 298.15, P1, q
h0 = gas.enthalpy_mole*10**-6

gas.TPX = 1E-14, P1, q

f = open('debug_enthalpy.txt', 'w')
f.write('{:4s}\n\n'.format(spec))


h = np.zeros_like(T)
s = np.zeros_like(T)
for i in range(0,len(T)):
    gas.TP = T[i], P1
    h[i] = gas.enthalpy_mole*10**-6 - h0
    s[i] = gas.entropy_mole*10**-3

    f.write('{:<12.0f}{:>12.6f}{:>20.6f}\n'.format(T[i], h[i], s[i]))
    print('{:<12.2f}{:>12.6f}{:>20.6f}'.format(T[i], h[i], s[i]))