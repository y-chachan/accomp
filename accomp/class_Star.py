import numpy as np
from pkg_resources import resource_filename

from .utils import get_elem_ratio, get_mass_fraction, get_elem_metal_mass_fraction

class Star:
    def __init__(self, mass=1., ref_abund='Asplund2021', custom_abund=None):
        """        
        Creates an instance of Star object.

        Parameters
        ----------
        mass: stellar mass in solar mass
        ref_abund: Load an abundance dictionary from Asplund 2009, Asplund 2021, or Lodders 2025
        custom_abund: supply a customized stellar abundance file. Format must follow that of pre-set reference abundance files.
        """  

        self.mass = mass

        if custom_abund is None:
            assert(ref_abund is not None)
            if ref_abund == 'Asplund2009':
                ref_abund_file = resource_filename(__name__, "data/asplund_2009_solar_abundance.in")
            elif ref_abund == 'Asplund2021':
                ref_abund_file = resource_filename(__name__, "data/asplund_2021_solar_abundance.in")
            elif ref_abund == 'Lodders2025':
                ref_abund_file = resource_filename(__name__, "data/lodders_2025_solar_abundance.in")
            else:
                print('Specified composition is not currently implemented.')
                raise ValueError

        if custom_abund is not None:
            assert(ref_abund is None)
            ref_abund_file = custom_abund

        self.abund = np.genfromtxt(ref_abund_file, usecols=(1))
        self.species_name = np.genfromtxt(ref_abund_file, usecols=0, dtype=str)
        self.abundance_dict = {}

        for i, name in enumerate(self.species_name):
            self.abundance_dict[name] = self.abund[i]

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
        """returns a bare bones abundance dictionary with H and He. Useful for adding elements or creating custom abundance dictionaries."""
        return {'H':self.abundance_dict['H'], 'He':self.abundance_dict['He']}

    def get_stellar_normalized_abundance(self, abund_dict):
        """
        Function for converting log abundances to metallicity relative to stellar. Takes a dictionary of abundances as input.
        """
        metallicity_times_stellar = {}
        for sp in abund_dict.keys():
            if sp in ['H', 'He']:continue
            metallicity_times_stellar[sp] = 10**(abund_dict[sp] - self.abundance_dict[sp])
        return metallicity_times_stellar

    def convert_metallicity_to_abundance(self, metallicity_dict):
        """
        Function for converting metallicity into log abundance. Takes a dictionary of metallicity as input.        
        """
        abd_dict = self.get_basic_abund_dict()
        for sp in metallicity_dict.keys():
            abd_dict[sp] = self.abundance_dict[sp] + np.log10(metallicity_dict[sp])

        return abd_dict

