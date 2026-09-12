import numpy as np
import textwrap
from .utils import get_mass_fraction

class Planet:
    def __init__(self, star, measured_abund_dict, loc=None):
        """
        Creates an instance of Planet object.

        Parameters
        ----------
        star_object: input an instance of class_Star for stellar properties. The stellar composition assigned to this object must correspond to the values relative to which the abundance dictionary for the planet is provided.
        measured_abund_dict: must be a dictionary containing abundance of different elements with respect to H normalized to stellar values. The keys of the dictionary must correspond to the different elements that have been measured.
        loc: float specifying the location of a planet in au
        """

        self.star = star
        self.measured_abund = measured_abund_dict
        self.loc = loc

        self.species_name = self.measured_abund.keys()
        self.set_abundance_dict()
        self.mass_averaged_enrichment = self.get_mass_averaged_enrichment(self.measured_abund)
        self.mass_fraction = get_mass_fraction(self.abund_dict) #X, Y, Z

    def set_abundance_dict(self):
        """Set up an abundance dictionary with log n_H = 12 and additional elements that have been measured."""
        self.abund_dict = self.star.get_basic_abund_dict()

        for sp in self.species_name:
            self.abund_dict[sp] = self.star.abundance_dict[sp] + np.log10(self.measured_abund[sp]) 

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
        """Returns the mass fraction averaged enrichment of elements in the planet."""

        mass_averaged_enrichment = 0.
        for sp in enrichment_dict.keys():
            if sp in ['H', 'He']:continue
            mass_averaged_enrichment += enrichment_dict[sp] * self.star.Zi_over_Z[sp]
        return mass_averaged_enrichment

    def get_component_mass_fraction(self, component_enrichment_dict):
        """Returns the mass fraction contribution of a subset of metals accreted through a given source, e.g., local gas, solids, or enriched gas."""

        return self.get_mass_averaged_enrichment(component_enrichment_dict) / self.mass_averaged_enrichment * self.mass_fraction[2]


    @staticmethod
    def simulate_refractory_measurement(measured_abund, most_ref_element, refractory_list, solids_comp):
        """For a given composition of disk solids, take the element with largest fraction in solids amongst the measured ones and use it to calculate the enrichment we would expect for a fully refractory species. """

        simulated_measurements = {}
        for ref in refractory_list:
            if ref in measured_abund.keys():
                print('Element ' + ref + ' already amongst measured species.')
                continue
            simulated_measurements[ref] = (measured_abund[most_ref_element] - (1 - solids_comp[most_ref_element])) / solids_comp[most_ref_element]
        return simulated_measurements

    def get_metal_mass_breakdown(self, solids_comp, gas_comp, main_ref_elem):
        """Pass in a dictionary of the solid and gas composition and obtain the metal mass accreted via local solids and gas, and excess metals that are unnacounted for by local material.
    
        Parameters
        ----------   
        solids_comp: A dictionary specifying the fraction of each element in solid phase.
        gas_comp: A dictionary specifying the fraction of each element in gas phase.
        main_ref_elem: The element amongst the measured ones that will constitute the refractory reference. 
        """

        self.enrichment_from_solid = {}
        self.enrichment_from_local_gas = {}
        self.enrichment_from_excess_metal = {}
        
        #The following equations assume this.
        assert(solids_comp[main_ref_elem] == 1.), print(textwrap.fill('The reference refractory element must be present entirely in solids. If no such element exists amongst the measured ones, use method simulate_refractory_measurement to simulate measurements for a fully refractory species and add it to planet measurements'))

        for elem in self.species_name:
            if solids_comp[elem] == 1.:
                self.enrichment_from_solid[elem] = self.measured_abund[elem]
                self.enrichment_from_local_gas[elem] = 0.
            else:
                self.enrichment_from_solid[elem] = self.measured_abund[main_ref_elem] * solids_comp[elem]
                self.enrichment_from_local_gas[elem] = gas_comp[elem]
                self.enrichment_from_excess_metal[elem] = (self.measured_abund[elem] - 1) - (self.measured_abund[main_ref_elem] - 1) * solids_comp[elem] #np.maximum((self.measured_abund[elem] - 1) - (self.measured_abund[main_ref_elem] - 1) * solids_comp[elem], 0.)

            if np.any(self.enrichment_from_local_gas[elem] + self.enrichment_from_solid[elem] > self.measured_abund[elem]):
                print(textwrap.fill('For a fraction of samples, the predicted enrichment for element ' + elem + ' is greater than the measured enrichment due to the assumed split for this element between local solid and gas composition. This is expected because planet measurements can \'fall\' on both sides of the predicted composition. Use function constrain_elem_frac_in_solids to put empirical constraints on the split of elements between disk solids and gas.'))

        self.mass_fraction_from_solids = self.get_component_mass_fraction(self.enrichment_from_solid)
        self.mass_fraction_from_local_gas = self.get_component_mass_fraction(self.enrichment_from_local_gas)
        self.mass_fraction_from_excess_metal = self.get_component_mass_fraction(self.enrichment_from_excess_metal)

        #convert mass fraction to absolute mass
        assert(hasattr(self, 'metal_mass'))
        self.metal_mass_solids = self.mass_fraction_from_solids * self.mass
        self.metal_mass_local_gas = self.mass_fraction_from_local_gas * self.mass
        self.metal_mass_excess_metal = self.mass_fraction_from_excess_metal * self.mass

    def constrain_elem_frac_in_solids(self, main_ref_elem):
        """Takes the measured abundances to estimate the fraction of each element that must be present in the solids such that measured abundance is a result of accretion of local solids and gas.
    
        Parameters
        ----------   
        main_ref_elem: The element amongst the measured ones that will constitute the refractory reference. 
        """

        f_v_dict = {}

        for elem in self.species_name:
            f_v_dict[elem] = (self.measured_abund[elem] - 1.) / (self.measured_abund[main_ref_elem] - 1.)

        return f_v_dict