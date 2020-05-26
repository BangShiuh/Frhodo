# This file is part of Frhodo. Copyright © 2020, UChicago Argonne, LLC
# and licensed under BSD-3-Clause. See License.txt in the top-level 
# directory for license and copyright information.

import numpy as np
import cantera as ct
import warnings
from scipy.optimize import curve_fit, OptimizeWarning


def fit_rate_eqn(P, X, mech, coefNames, rxnIdx):
    rxn = mech.gas.reaction(rxnIdx)
    def inner(invT, *coeffs):
        if type(rxn) is ct.ElementaryReaction or type(rxn) is ct.ThreeBodyReaction: # if 2 coeffs for Arrhenius
            if len(coeffs) == 2:                                                    # assume n = 0
                coeffs = np.append(coeffs, 0)

        for coefName, coefVal in zip(coefNames, coeffs):   # updates reaction mechanism for specific reaction
            mech.coeffs[rxnIdx][coefName] = coefVal
        mech.modify_reactions(mech.coeffs, rxnNums=rxnIdx)
        
        rate = []
        temperatures = np.divide(10000, invT)
        for n, T in enumerate(temperatures):  # TODO FOR PRESSURE DEPENDENT REACTIONS ZIP PRESSURE
            mech.set_TPX(T, P[n], X[n])
            rate.append(mech.gas.forward_rate_constants[rxnIdx])
        return np.log10(rate)
    return inner

def jacobian(func, x, h=np.finfo(float).eps): # central finite difference
    def OoM(x):
        x = np.copy(x)
        if not isinstance(x, np.ndarray):
            x = np.array(x)
        x[x==0] = 1                       # if zero, make OoM 0
        return np.floor(np.log10(np.abs(x)))
      
    len_x = len(x)
    h = np.ones((len_x, 1))*h
    h = (h.T * np.power(10, OoM(x))).T                  # scale h to OoM of variable
    x = np.tile(x, (2*len_x,1))
    x[1::2] = x[1::2] + (np.eye(len_x, len_x).T * h).T  # add h on odd rows
    x[::2] = x[::2] - (np.eye(len_x, len_x).T * h).T    # subtract h on even rows
    df = np.empty((len_x, 1))
    for i in range(len_x):
        df[i] = (func(x[2*i+1])[0] - func(x[2*i])[0])/(2*h[i])

    return df.T[0]

def scale_fcn(coefVals, coefNames, rxn, dir='forward'):
    coefVals = np.copy(coefVals)
    if type(rxn) is ct.ElementaryReaction or type(rxn) is ct.ThreeBodyReaction: # TODO: UPDATE THIS FOR OTHER TYPES OF EXPRESSIONS
        for n, coefVal in enumerate(coefVals):   # updates reaction mechanism for specific reaction
            if coefNames[n] == 'pre_exponential_factor':
                if dir == 'forward':
                    coefVals[n] = 10**coefVal
                else:
                    coefVals[n] = np.log10(coefVal)
    return coefVals

def fit_coeffs(rates, T, P, X, coefNames, rxnIdx, mech):
    if len(coefNames) == 0: return # if not coefs being optimized in rxn, return 
    
    min_neg_system_value = np.finfo(float).min*(1E-20) # Don't push the limits too hard
    min_pos_system_value = np.finfo(float).eps*(1.1)
    max_pos_system_value = np.finfo(float).max*(1E-20)
    
    rxn = mech.gas.reaction(rxnIdx)
    x0 = []
    lower_bnd = []
    upper_bnd = []
    for n, coefName in enumerate(coefNames):
        # if coef['rxnIdx'] != rxnIdx: continue   # ignore reaction not specified
        x0.append(mech.coeffs_bnds[rxnIdx][coefName]['resetVal'])
        if np.isnan(mech.coeffs_bnds[rxnIdx][coefName]['limits']).any():
            if coefName == 'pre_exponential_factor':
                lower_bnd.append(min_pos_system_value)             # A should be positive
            elif coefName == 'activation_energy' and x0[n] > 0:
                lower_bnd.append(0)                                # Ea shouldn't change sign
            else:
                lower_bnd.append(min_neg_system_value)
            
            if coefName == 'activation_energy' and x0[n] < 0:   # Ea shouldn't change sign
                upper_bnd.append(0)
            else:
                upper_bnd.append(max_pos_system_value)
        else:
            lower_bnd.append(mech.coeffs_bnds[rxnIdx][coefName]['limits'][0])
            upper_bnd.append(mech.coeffs_bnds[rxnIdx][coefName]['limits'][1])
        
    x0s = scale_fcn(x0, coefNames, rxn, dir='inverse')
        
    invT = np.divide(10000, T)
    logk = np.log10(rates)
    
    if not isinstance(X, (list, np.ndarray)):   # if only a single composition is given, duplicate
        X = [X]*len(invT)
    
    eqn = lambda invT, *x: fit_rate_eqn(P, X, mech, coefNames, rxnIdx)(invT, *scale_fcn(x, coefNames, rxn))
    s = np.abs(jacobian(lambda x: eqn([np.mean(invT)], *x), x0s, 1E-9))
    s[s==0] = 1E-9  # TODO: MAKE THIS BETTER running into problem when Ea is zero, this is a janky workaround
    # s /= np.max(s)  # to prevent overflow if s_i is > 1 and unbounded
    scaled_eqn = lambda invT, *x: eqn(invT, *(x/s + x0s))
    lower_bnd = (scale_fcn(lower_bnd, coefNames, rxn, dir='inverse') - x0s)*s
    upper_bnd = (scale_fcn(upper_bnd, coefNames, rxn, dir='inverse') - x0s)*s
    p0 = np.zeros_like(x0s)
      
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', OptimizeWarning)
        try:
            popt, _ = curve_fit(scaled_eqn, invT, logk, p0=p0, bounds=[lower_bnd, upper_bnd],
                            method='dogbox')
        except:
            return
    
    coeffs = scale_fcn(popt/s + x0s, coefNames, rxn)
    # print(coeffs)
    # print('')
    
    # T_test = np.linspace(np.min(T)*0.999, np.max(T)*1.001, 50)
    # rate_test = []
    # for T_t in T_test:
        # mech.set_TPX(T_t, P[0], X[0])
        # rate_test.append(mech.gas.forward_rate_constants[rxnIdx])
    
    # import matplotlib.pyplot as plt     
    # plt.plot(np.divide(10000, T), np.log10(rates), 'o')
    # plt.plot(np.divide(10000, T_test), np.log10(rate_test))
    # plt.show()
    
    return coeffs

if __name__ == '__main__':
    import mech_fcns
    mech = mech_fcns.Chemical_Mechanism()
    try:
        mech.load_mechanism({'mech': 'generated_mech.yaml'})
    except:
        print("generated_mech.yaml doesn't exist in the current directory")
        quit()
    
    # rates = [1529339.05689338, 1548270.86688399, 1567437.0352583]
    rates = [1529339.05689338, 1548270.86688399, 1567437.0352583]*np.array([1, 1.001, 1])
    T = [2387.10188629, 2389.48898818, 2391.88086905]
    P = [16136.20900077, 16136.20900077, 16136.20900077]
    X = {'Kr': 0.99, 'C8H8': 0.01}
    
    coefNames = ['activation_energy', 'pre_exponential_factor', 'temperature_exponent']
    rxnIdx = 0
    mech = self.parent.mech
    fit_coeffs(rates, T, P, X, coefNames, rxnIdx, mech)
    print(np.array([2.4442928e+08, 3.4120000e+11, 0.0000000e+00]))