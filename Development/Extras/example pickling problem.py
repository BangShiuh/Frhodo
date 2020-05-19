import cantera as ct
import multiprocessing as mp

def parallel_func(x):
    IC, gas = x
    # do things with cantera
    return IC['T']  # return output value from cantera

gas = ct.Solution('gri30.xml')
gas.X = {'Ar': 0.99, 'O2': 0.003, 'H2': 0.007}

IC_list = [{'T': 1200, 'P': 8000}, {'T': 1500, 'P': 8000}]

output = []
pool = mp.Pool(1)
output = pool.map_async(parallel_func, [(IC, gas) for IC in IC_list]).get()
pool.close()