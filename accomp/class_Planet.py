import numpy as np
from .utils import get_mass_fraction

class Planet:
    def __init__(self, star, measured_abund_dict):
        """
        Creates an instance of Planet object.

        Parameters
        ----------
        star_object: input an instance of class_Star for stellar properties. The stellar composition assigned to this object must correspond to the values relative to which the abundance dictionary for the planet is provided.
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
        Parameters
        ----------   
        mass: Either a numpy array or a float. If numpy array, must be the same length as chain containing planet abundance.
        mass_err: If none, just use mass. If a float, mass must be a float too, will generate a random normal distribution sample with length the same as planet abundance chains. 
        """

        if mass_err is not None:
            assert(type(mass) == float)
            self.mass = np.random.normal(mass, mass_err, len(self.mass_fraction[2]))
        else:
            self.mass = mass

    def get_metal_mass(self):
        """Returns the absolute amount of metals in the planet. Same unit as planet mass. Set planet mass before calling this function"""

        assert(hasattr(self, 'mass'))
        self.metal_mass = self.mass_fraction[2] * self.mass

    def get_mass_averaged_enrichment(self, enrichment_dict):
        """Returns the absolute amount of metals in the planet. Same unit as planet mass. Set planet mass before calling this function"""

        mass_averaged_enrichment = 0.
        for sp in enrichment_dict.keys():
            if sp in ['H', 'He']:continue
            mass_averaged_enrichment += enrichment_dict[sp] * self.star.Zi_over_Z[sp]
        return mass_averaged_enrichment

    def get_component_mass_fraction(self, component_enrichment_dict):
        """Returns the mass fraction contribution of a subset of metals accreted through a given source, e.g., local gas, solids, or enriched gas."""

        return self.get_mass_averaged_enrichment(component_enrichment_dict) / self.mass_averaged_enrichment * self.mass_fraction[2]

    def get_metal_mass_breakdown(self, solids_comp, gas_comp, main_ref_elem):
        """Pass in a dictionary of the solid and gas composition and obtain the metal mass accreted via local solids and gas, and excess metals that are unnacounted for by local material.
    
        Parameters
        ----------   
        solids_comp: A dictionary specifying the fraction of each element in solid phase.
        gas_comp: A dictionary specifying the fraction of each element in gas phase.
        main_ref_elem: The element amongst the measured ones that will constitute the refractory reference. 
        """

        self.enrichment_from_solid = self.star.get_basic_abund_dict()
        self.enrichment_from_local_gas = self.star.get_basic_abund_dict()
        self.enrichment_from_excess_metal = self.star.get_basic_abund_dict()

        #The following equations assume this.
        assert(solids_comp[main_ref_elem] == 1.), print('The reference refractory element must be present entirely in solids.')

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