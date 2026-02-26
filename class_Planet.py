import numpy as np
from utils import get_mass_fraction

class Planet:
    def __init__(self, star, measured_abund_dict):
        """
        star: object of class Star and contains the stellar composition relative to which the abundance dictionary for the planet is provided

        measured_abund_dict: must be a dictionary containing abundance of different elements with respect to H normalized to stellar values. The keys of the dictionary must correspond to the different elements that have been measured.
        """
        self.star = star
        self.measured_abund = measured_abund_dict
        self.species_name = self.measured_abund.keys()
        self.abund_dict = self.star.get_basic_abund_dict()

        for sp in self.species_name:
            self.abund_dict[sp] = star.abundance_dict[sp] + np.log10(self.measured_abund[sp]) 

        self.mass_averaged_enrichment = self.get_mass_averaged_enrichment(self.measured_abund)
        self.mass_fraction = get_mass_fraction(self.abund_dict) #X, Y, Z


    def set_planet_mass(self, mass, mass_err=None):
        """        
        :param mass: Either a numpy array or a float. If numpy array, must be the same length as chain containing planet abundance.
        :param mass_err: If none, just use mass. If a float, mass must be a float too, will generate a random normal distribution sample with length the same as planet abundance chains. 
        """
        if mass_err is not None:
            assert(type(mass) == float)
            self.mass = np.random.normal(mass, mass_err, len(self.mass_fraction[2]))
        else:
            self.mass = mass

    def get_metal_mass(self):
        """Set planet mass before calling this function"""
        assert(hasattr(self, 'mass'))
        self.metal_mass = self.mass_fraction[2] * self.mass

    def get_mass_averaged_enrichment(self, enrichment_dict):
        mass_averaged_enrichment = 0.
        for sp in enrichment_dict.keys():
            if sp in ['H', 'He']:continue
            mass_averaged_enrichment += enrichment_dict[sp] * self.star.Zi_over_Z[sp]
        return mass_averaged_enrichment

    def get_component_mass_fraction(self, component_enrichment_dict):
        return self.get_mass_averaged_enrichment(component_enrichment_dict) / self.mass_averaged_enrichment * self.mass_fraction[2]

    def get_metal_mass_breakdown(self, solids_comp, gas_comp, main_ref_elem):
        """Pass in an indication of the solid mass fraction of each element and obtain the mass accreted via solids and local gas"""
        self.enrichment_from_solid = self.star.get_basic_abund_dict()
        self.enrichment_from_local_gas = self.star.get_basic_abund_dict()
        self.enrichment_from_excess_metal = self.star.get_basic_abund_dict()

        for elem in self.species_name:
            if solids_comp[elem] == 1.:
                self.enrichment_from_solid[elem] = self.measured_abund[elem]
            else:
                self.enrichment_from_solid[elem] = self.measured_abund[main_ref_elem] * solids_comp[elem]
                self.enrichment_from_local_gas[elem] = gas_comp[elem]
                self.enrichment_from_excess_metal[elem] = np.maximum((self.measured_abund[elem] - 1) - (self.measured_abund[main_ref_elem] - 1) * solids_comp[elem], 0.)

        self.mass_fraction_from_solids = self.get_component_mass_fraction(self.enrichment_from_solid)
        self.mass_fraction_from_local_gas = self.get_component_mass_fraction(self.enrichment_from_local_gas)
        self.mass_fraction_from_excess_metal = self.get_component_mass_fraction(self.enrichment_from_excess_metal)

        #convert mass fraction to absolute mass
        assert(hasattr(self, 'metal_mass'))
        self.metal_mass_solids = self.mass_fraction_from_solids * self.mass
        self.metal_mass_local_gas = self.mass_fraction_from_local_gas * self.mass
        self.metal_mass_excess_metal = self.mass_fraction_from_excess_metal * self.mass


    #TODO: write module to compare planetary measurements with expectations