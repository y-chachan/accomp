import numpy as np
import astropy.constants as const
from .class_Molecule import Molecule, EnrichedMolecule, MoleculeDict
from scipy.interpolate import interp1d

mu = 2.3 #mean molecular weight of disk gas

class Disk:
    def __init__(self, star_object, n_r=30, r_min=0.1, r_max=200., alpha=1e-3, d2g=0.01):
        """        
        Creates an instance of Disk object.

        Parameters
        ----------
        star_object: input an instance of class_Star for stellar properties
        n_r: number of radial grid points
        r_min: inner edge of the disk in au
        r_max: outer edge of the disk in au
        alpha: dimensionless viscosity parameter
        d2g: dust to gas mass ratio of the disk
        """  

        self.star = star_object
        self.mstar = self.star.mass * const.M_sun.cgs.value
        self.semimajor = np.logspace(np.log10(r_min), np.log10(r_max), n_r) * const.au.cgs.value
        self.Omega_K = np.power(const.G.cgs.value * self.mstar / self.semimajor**3, 1/2)
        self.per = 2 * np.pi / self.Omega_K
        self.d2g = 0.01
        self.alpha = alpha

    def set_temp_struct(self, M_dot_solar=1e-8, dlnMdotdlnM=2., dlnLdlnM=1.5):
        """        
        Sets the disk's temperature profile.

        Parameters
        ----------
        M_dot_solar: accretion rate onto star in units of solar masses/year
        dlnMdotdlnM: power law exponent relating stellar accretion rate to stellar mass
        dlnLdlnM: power law exponent relating stellar luminosity to stellar mass
        """    

        self.Mdot_mass_exponent = dlnMdotdlnM
        self.M_dot = M_dot_solar * (self.mstar/const.M_sun.cgs.value)**dlnMdotdlnM

        alpha_fiducial=1e-3
        M_dot_fiducial=1e-8

        M_dot_scaled = (self.M_dot / M_dot_fiducial) #relative to fiducial
        L_star_scaled = (self.mstar/const.M_sun.cgs.value)**dlnLdlnM #relative to fiducial

        T_vis = 200. * (self.mstar/const.M_sun.cgs.value)**0.3 * (self.alpha/alpha_fiducial)**-0.2 * (M_dot_scaled)**0.4 * (self.semimajor/const.au.cgs.value)**-0.9
        T_irr = 150. * (L_star_scaled)**(2/7) * (self.mstar/const.M_sun.cgs.value)**(-1/7) * (self.semimajor/const.au.cgs.value)**(-3/7)
        self.T_disk = np.minimum(2000, np.maximum(T_vis, T_irr)) #assuming temperature does not exceed 2000

        #disk temperature allows us to calculate sound speed, scale height, and gamma=dlnP/dlnr
        self.c_s = np.sqrt(const.k_B.cgs.value * self.T_disk / mu / const.m_p.cgs.value)
        self.H_g = self.c_s / self.Omega_K
        self.gamma = np.gradient(np.log(self.Omega_K**2 / self.c_s), np.log(self.semimajor))

    def get_pebble_Miso(self):
        """calculate the pebble isolation mass based on Bitsch et al. 2018. Requires the disk temperature profile to be set in advance."""

        assert(hasattr(self, 'T_disk')), "Need to set disk temperature profile first"

        self.pebble_Miso = 25. * (self.mstar/const.M_sun.cgs.value) * (self.H_g / self.semimajor / 0.05)**3 * (0.34 * (-3 / np.log10(self.alpha))**4 + 0.66) * (1 - (self.gamma + 2.5)/6)

    def get_temperature_location(self, temp):
        """obtain the semimajor axis corresponding to a specifed temperature. Useful for calculating snowline location"""

        return 10**interp1d(np.log10(self.T_disk), np.log10(self.semimajor/const.au.cgs.value))(np.log10(temp))

    def get_location_temperature(self, distance):
        """obtain the temperature corresponding to a specifed semimajor axis in au."""

        return 10**interp1d(np.log10(self.semimajor/const.au.cgs.value), np.log10(self.T_disk))(np.log10(distance))

    def get_condensation_locations(self, mol):
        """obtain the location of condensation fronts(s) of a Molecule or molecules in MoleculeDict in au."""

        if isinstance(mol, Molecule):
            return self.get_temperature_location(mol.condensation_T)

        if isinstance(mol, MoleculeDict):
            mol_Tcond = mol.get_mol_Tcond()
            Tcond_loc = {}
            for m in mol_Tcond.keys():
                Tcond_loc[m] = self.get_temperature_location(mol_Tcond[m])

            return Tcond_loc

    def get_solid_gas_composition(self, mol_dict, temp=None, location=None):
        """        
        Return dictionaries of solid and gas composition at a specified temperature or location. Only one of them should be specified.

        Parameters
        ----------
        mol_dict: an instance of MoleculeDict with set composition
        temp: temperature at which the solid and gas composition is desired
        location: location at which the solid and gas composition is desired
        """

        assert(np.logical_xor(temp is not None, location is not None))

        if temp is None:
            temp = self.get_location_temperature(location)

        solids_dict, gas_dict = MoleculeDict(self.star), MoleculeDict(self.star)
        for m, mol in mol_dict.molecule_dict.items():
            if mol.condensation_T > temp:
                solids_dict.molecule_dict[m] = mol
            else:
                gas_dict.molecule_dict[m] = mol

        solids = dict(solids_dict.summed_abundances)
        gas = dict(gas_dict.summed_abundances)

        #this is done to ensure that all species are present in both solid and gas phase dictionaries
        keys = solids.keys() | gas.keys()
        solids = {k: solids.get(k, 0.) for k in keys}
        gas = {k: gas.get(k, 0.) for k in keys}

        return solids, gas
    
    def assemble_planet(self, ref_to_H, solid_comp, gas_comp):
        """
        Return the enrichment and abundance of all the elements for a given ref_to_H

        Parameters
        ----------
        ref_to_H: refractory to H ratio, equivalent to solid-to-gas accretion rate
        solid_comp: a dictionary containing the fraction of each element in solid phase 
        gas_comp: a dictionary containing the fraction of each element in gas phase  
        """

        enrichment_dict = {}
        for sp in solid_comp.keys():
            enrichment_dict[sp] = solid_comp[sp] * ref_to_H + gas_comp[sp]

        return enrichment_dict
