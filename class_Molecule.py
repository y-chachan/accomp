import numpy as np
import re
from collections import defaultdict, Counter
import utils

class Molecule():
    def __init__(self, formula, T_cond):
        self.formula = formula
        self.elem_dict = self.parse_formula(formula)
        self.mol_weight = self.calculate_molecular_weight()
        self.condensation_T = T_cond

    def parse_formula(self, formula: str) -> dict[str, int]:
        pattern = r"([A-Z][a-z]*)(\d*)"
        counts = defaultdict(int)

        for element, number in re.findall(pattern, formula):
            counts[element] += int(number) if number else 1

        return dict(counts)

    def calculate_molecular_weight(self):
        mol_weight = 0.

        for sp in self.elem_dict.keys():
            mol_weight += utils.atomic_weight_dict[sp] * self.elem_dict[sp]

        return mol_weight

    def set_elem_fraction(self, star, elem, elem_frac, enriched_molecule=False):
        """
        Given an element and the fraction of that element contained in this molecule, calculate the fraction of other elements contained in this molecule.
        
        star (stellar object): need to specify stellar composition because stellar ratios are needed to calculate the fraction for other elements
        elem (string): the element for which the fraction is specified
        elem_frac (float): the fraction of the given element contained in this molecule
        """
        self.star = star
        self.elem_frac_dict = {}
        self.elem_frac_dict[elem] = elem_frac

        for sp in self.elem_dict:
            if sp in [elem, 'H', 'He']:
                continue
            else:
                self.elem_frac_dict[sp] = elem_frac * self.elem_dict[sp] / self.elem_dict[elem] * utils.get_elem_ratio(star.abundance_dict, elem, sp)

        if enriched_molecule is False:
            if np.any(np.fromiter(self.elem_frac_dict.values(), dtype=float) > 1.):
                print(self.elem_frac_dict)
                print('Fraction of element in this molecule cannot exceed 1')
                raise Exception


class EnrichedMolecule(Molecule):
    def __init__(self, formula, T_cond):
        super().__init__(formula, T_cond)

    #TODO:write a function to specify spatially variable enrichment

    #is this the best way to do this?
    def set_elem_fraction(self, star, elem, enrichment):
        super().set_elem_fraction(star, elem, enrichment, enriched_molecule=True)


class MoleculeDict():
    def __init__(self, star):
        self.molecule_dict = {}
        self.star = star
        #self.set_default_composition()

    @property
    def summed_abundances(self):
        sum = Counter()

        enriched_molecules = False

        for mol in self.molecule_dict.values():
            sum += Counter(mol.elem_frac_dict)
            if isinstance(mol, EnrichedMolecule):
                enriched_molecules = True

        if any(s > 1. for s in sum.values()) & (enriched_molecules == False):
            print(sum)
            raise ValueError('Element fractions sum up to > 1. Analyze output.')
        elif any(s < 0. for s in sum.values()):
            print(sum)
            raise ValueError('Element fractions sum up to < 0. Analyze output.')
        else:
            return sum
            
    @property
    def elem_dict(self):
        out = defaultdict(dict)
        for outer_k, inner in self.molecule_dict.items():
            for inner_k, value in inner.elem_frac_dict.items():
                out[inner_k][outer_k] = value
        return dict(out)

    def get_mol_Tcond(self):
        Tcond_dict = {}
        for k, mol in self.molecule_dict.items():
            Tcond_dict[k] = mol.condensation_T

        return dict(sorted(Tcond_dict.items(), key=lambda x: x[1], reverse=True))

    def set_default_composition(self):
        """set up molecules for all the key species"""
        #C-bearing species
        self.molecule_dict['CO'] = Molecule('CO', 30.)
        self.molecule_dict['CO'].set_elem_fraction(self.star, 'C', 0.4)

        self.molecule_dict['CO2'] = Molecule('CO2', 70.)
        self.molecule_dict['CO2'].set_elem_fraction(self.star, 'C', 0.1)

        #N-bearing species
        self.molecule_dict['N2'] = Molecule('N2', 25.)
        self.molecule_dict['N2'].set_elem_fraction(self.star, 'N', 0.6)

        self.molecule_dict['NH4SH'] = Molecule('NH4SH', 300.) #check temperature
        self.molecule_dict['NH4SH'].set_elem_fraction(self.star, 'N', 0.15)

        self.molecule_dict['NH3'] = Molecule('NH3', 70.) #check temperature
        self.molecule_dict['NH3'].set_elem_fraction(self.star, 'N', 0.25)

        #Nobel gases
        for k in ['Ne', 'Ar']:
            self.molecule_dict[k] = Molecule(k, 25.)
            self.molecule_dict[k].set_elem_fraction(self.star, k, 1.)

        #O-bearing species
        self.molecule_dict['H2O'] = Molecule('H2O', 170.)
        self.molecule_dict['H2O'].set_elem_fraction(self.star, 'O', 0.45)

        #S-bearing species
        self.molecule_dict['H2S'] = Molecule('H2S', 90.)
        self.molecule_dict['H2S'].set_elem_fraction(self.star, 'S', 0.1)

        #refractories
        self.molecule_dict['Na'] = Molecule('Na', 1100.)
        self.molecule_dict['Na'].set_elem_fraction(self.star, 'Na', 1.)

        for k in ['Si', 'Mg', 'Fe', 'Ni']:
            self.molecule_dict[k] = Molecule(k, 1600.)
            self.molecule_dict[k].set_elem_fraction(self.star, k, 1.)

        for k in ['Al', 'Ca']:
            self.molecule_dict[k] = Molecule(k, 2000.)
            self.molecule_dict[k].set_elem_fraction(self.star, k, 1.)

        #catch remaining fraction - put it in unknown refractory species
        self.catch_remaining_fraction()


    def catch_remaining_fraction(self):
        summed_abundances = self.summed_abundances
        for sp, frac in summed_abundances.items():
            if (frac < 1.):
                print('Putting remaining element fraction for ' + sp + ' in refractory')
                if sp == 'O':
                    self.molecule_dict['O_refractory'] = Molecule('O', 1600.)
                    self.molecule_dict['O_refractory'].set_elem_fraction(self.star, 'O', 1. - summed_abundances['O'])
                elif sp == 'C':
                    self.molecule_dict['C_refactory'] = Molecule('C', 300.)
                    self.molecule_dict['C_refactory'].set_elem_fraction(self.star, 'C', 1. - summed_abundances['C'])
                elif sp == 'S':
                    #Refractory S in some form other than ammonium-sulfate salts
                    self.molecule_dict['S_refactory'] = Molecule('S', 300.)
                    self.molecule_dict['S_refactory'].set_elem_fraction(self.star, 'S', 1. - summed_abundances['S'])
                elif sp == 'N':
                    #Refractory N in some form other than ammonium-sulfate salts
                    self.molecule_dict['N_refactory'] = Molecule('N', 300.)
                    self.molecule_dict['N_refactory'].set_elem_fraction(self.star, 'N', 1. - summed_abundances['N'])
                else:
                    print('Species' + sp + ' not implemented yet.' )

