import numpy as np
from scipy import constants as const
k_ev = const.k / const.eV
T_def = 298.15

def pka_to_delta_g(pka, n_protons=1, T=T_def):
    return k_ev * T * np.log(10) * pka * n_protons

def delta_g_to_pka(delta_g, n_protons=1, T=T_def):
    return delta_g / (k_ev * T * np.log(10) * n_protons)

data_dict = {
    "NO3H": -1.3,
    "NO2H": 3.29,
    "H3O+": 0.0,
    "H2O": 14.0,
    "NH4+": 9.24,
}