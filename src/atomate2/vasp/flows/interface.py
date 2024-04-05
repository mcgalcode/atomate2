
from ..jobs.core import RelaxSetGenerator, RelaxMaker
from pymatgen.core import Structure
from pymatgen.analysis.elasticity import Strain

from jobflow import Flow

class InterfaceEnergyMaker():


    def make(self, struct1, struct2):

        #optimize thickness by iterating up number of layers
        # optimize kpoints if insulating?


        relax_jobs = [self.get_relax_job(s) for s in [struct1, struct2]]

        interface_struct, strains = self.calculate_interface_struct(struct1, struct2)

        interface_energy_job = self.get_interface_relax_job(interface_struct)
        #HOI: what is the interface relax job?
        strained_bulk_surface_jobs = [self.get_strained_bulk_job(s, st) for s, st in zip([struct1, struct2], strains)]

        return Flow([*relax_jobs, interface_energy_job, *strained_bulk_surface_jobs])
    
    def get_strained_bulk_job(self, structure: Structure, strain: Strain, termination):
    
        strained_structure = self.apply_strain(structure, strain)

        strained_surface = self.get_surface(strained_structure, termination)
    
        input_set_generator = RelaxSetGenerator(user_incar_settings={"ISIF": 2})
        return RelaxMaker()
    
    def apply_strain(self, struct, strain):
        pass

    def get_surface(self, struct, termination):
        pass