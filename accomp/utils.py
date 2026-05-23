import numpy as np
from pkg_resources import resource_filename

species_name = np.loadtxt(resource_filename(__name__, "data/species_info.txt"), usecols=1, dtype=str)
conv = lambda x: float(x.decode().split('(')[0]) #this is to remove the uncertainty in parenthesis before importing
atomic_weight = np.loadtxt(resource_filename(__name__, "data/species_info.txt"), usecols=3, converters=conv)
atomic_weight_dict = {}
for i, name in enumerate(species_name):
    atomic_weight_dict[name] = atomic_weight[i]

#TODO: we also need something to indicate whether a species is volatile or refractory


def get_elem_ratio(abund_dict, e1, e2):
    return np.power(10, abund_dict[e1] - abund_dict[e2]) 
   
   
def get_mass_fraction(abundance_dict):
    X = 10**abundance_dict['H'] * atomic_weight_dict['H']
    Y = 10**abundance_dict['He'] * atomic_weight_dict['He']
    Z = 0.
    
    for sp in abundance_dict.keys():        
        if sp in ['H', 'He']:
            continue
        else:
            Z += 10**abundance_dict[sp] * atomic_weight_dict[sp]

    normalize_factor = X + Y + Z
    return X / normalize_factor, Y / normalize_factor, Z / normalize_factor


def get_elem_metal_mass_fraction(abundance_dict):
    X, Y, Z = get_mass_fraction(abundance_dict)

    Zi_over_X = {}
    Zi_over_Z = {}

    for sp in abundance_dict.keys():
        if sp in ['H', 'He']:
            continue
        else:
            Zi_over_X[sp] = 10**abundance_dict[sp] * atomic_weight_dict[sp] / (10**abundance_dict['H'] * atomic_weight_dict['H'])
            Zi_over_Z[sp] = Zi_over_X[sp] / (Z/X)

    return Zi_over_X, Zi_over_Z