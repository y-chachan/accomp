import numpy as np
from utils import get_elem_ratio, get_mass_fraction, get_elem_metal_mass_fraction

class Star:
    def __init__(self, mass=1., ref_abund='asplund_2021_solar_abundance.in'):
        self.mass = mass
        self.abund = np.genfromtxt(ref_abund, usecols=(1,2))
        self.species_name = np.genfromtxt(ref_abund, usecols=0, dtype=str)
        self.abundance_dict = {}

        for i, name in enumerate(self.species_name):
            self.abundance_dict[name] = self.abund[i,0]

        self.CtoO = get_elem_ratio(self.abundance_dict, 'C', 'O')
        self.OtoSi = get_elem_ratio(self.abundance_dict, 'O', 'Si')
        self.CtoSi = get_elem_ratio(self.abundance_dict, 'C', 'Si')
        self.num_species = len(self.abundance_dict.keys())
        self.X, self.Y, self.Z = get_mass_fraction(self.abundance_dict)

        self.Z_over_X = self.Z / self.X
        self.X_prime = self.X / (1. - self.Z)
        self.Y_prime = 1. - self.X_prime
        self.Zi_over_Z, self.Zi_over_Z = get_elem_metal_mass_fraction(self.abundance_dict)

    def get_basic_abund_dict(self):
        return {'H':self.abundance_dict['H'], 'He':self.abundance_dict['He']}

    def get_stellar_normalized_abundance(self, abund_dict):
        """
        Function for converting log abundances to metallicity relative to stellar        
        :param abund_dict: Description
        """
        metallicity_times_stellar = {}
        for sp in abund_dict.keys():
            if sp in ['H', 'He']:continue
            metallicity_times_stellar[sp] = 10**(abund_dict[sp] - self.abundance_dict[sp])
        return metallicity_times_stellar

    def convert_metallicity_to_abundance(self, metallicity_dict):
        """
        Function for converting metallicity into log abundance        
        :param metallicity_dict: Description
        """
        abd_dict = self.get_basic_abund_dict()
        for sp in metallicity_dict.keys():
            abd_dict[sp] = self.abundance_dict[sp] + np.log10(metallicity_dict[sp])

        return abd_dict

