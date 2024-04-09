import numpy as np

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from jobflow import Flow, Maker
from atomate2.vasp.jobs.core import RelaxMaker
from atomate2.vasp.sets.core import RelaxSetGenerator
from pymatgen.io.vasp import Kpoints
from pymatgen.core import Structure
from pymatgen.analysis.elasticity import Strain
from pymatgen.analysis.interfaces.substrate_analyzer import SubstrateAnalyzer
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder

if TYPE_CHECKING:
    from atomate2.vasp.jobs.core import BaseVaspMaker

class InterfaceEnergyMaker(Maker):
    """
    Maker to calculate the interface energy between two bulk structures.

    Parameters
    ----------
    name : str
        Name of the flows produced by this maker.
    relax_maker : .BaseVaspMaker
        Maker to produce the relaxation jobs for the structures.
    """

    name: str = "interface energy"
    relax_maker: BaseVaspMaker = field(default_factory=RelaxMaker)

    def get_relax_job(self, structure: Structure):
        input_set_generator = RelaxSetGenerator(user_incar_settings={"ISIF": 2})
        return self.relax_maker.make(structure, input_set_generator, user_kpoint_settings=Kpoints())
    
    #def get_surface(self, structure: Structure, termination: tuple):
    #    pass
    
    #def get_strained_surface(self, structure: Structure, strain: Strain, termination=None):
    #    strained_structure = structure.apply_strain(np.diag(strain))
        #strained_surface = self.get_surface(strained_structure, termination)
    #    return strained_structure
    
    def calculate_interface_structure(self, structure_1: Structure, structure_2: Structure):
        """
        Logic: 
        Two possible cases: structure_1 is the film and structure_2 is the substrate, or vice versa; film is the one that is strained.
        1. Calculate the substrate-analyzer match for both cases.
        2. Build the interface structure using the coherent interface builder, one interface for each termination; take only the first 5(?) terminations.
        3. Return the interface structure and the strains that were applied to the bulk structures. 
        Note: Have to account for the fact that interfacial energy is E_interface - [E(strained structure_1) + E(unstrained structure_2)] and vice versa.
        Also, average calculated energies for all termination sequences.
        """

        sub_analyzer = SubstrateAnalyzer()
        match_1 = list(sub_analyzer.calculate(structure_1, structure_2, lowest=True))
        match_2 = list(sub_analyzer.calculate(structure_2, structure_1, lowest=True))
        interfaces = []

        cib_1 = CoherentInterfaceBuilder(film_structure=structure_1,
                                substrate_structure=structure_2,
                                film_miller=match_1[0].film_miller,
                                substrate_miller=match_1[0].substrate_miller)
        for index, termination in enumerate(cib_1.terminations):
            if index == 5:
                break
            interfaces.append(list(cib_1.get_interfaces(termination=termination, gap=2.0, vacuum_over_film=0, film_thickness=3, substrate_thickness=1)))

        cib_2 = CoherentInterfaceBuilder(film_structure=structure_2,
                                substrate_structure=structure_1,
                                film_miller=match_2[0].film_miller,
                                substrate_miller=match_2[0].substrate_miller)
        for index, termination in enumerate(cib_2.terminations):
            if index == 5:
                break
            interfaces.append(list(cib_2.get_interfaces(termination=termination, gap=2.0, vacuum_over_film=0, film_thickness=3, substrate_thickness=1)))
        
        strained_structure_1 = structure_1.apply_strain(np.diag(match_1[0].strain)) #get_strained_surface(structure_1, match_1[0].strain, cib_1.terminations[0])
        strained_structure_2 = structure_1.apply_strain(np.diag(match_1[0].strain)) #get_strained_surface(structure_2, match_2[0].strain, cib_2.terminations[0])

        return interfaces, [strained_structure_1, strained_structure_2]
    
    def make(self, structure_1: Structure, structure_2: Structure):
        #TODO: 
        #1. optimize thickness by iterating up number of layers when building interface
        #2. optimize kpoints if insulating?
        #3. output job to calculate interface energy?
        #4. are surface jobs needed? Accordingly modify the get_strained_surface method.

        relax_jobs = [self.get_relax_job(s) for s in [structure_1, structure_2]]

        interfaces, strained_structures = self.calculate_interface_structure(structure_1, structure_2)

        interface_energy_jobs = [self.get_relax_job(interface) for interface in interfaces]

        strained_bulk_surface_jobs = [self.get_relax_job(strained_structure) for strained_structure in strained_structures]

        #output_job = self.get_output_job()?????

        return Flow([*relax_jobs, *interface_energy_jobs, *strained_bulk_surface_jobs])