import numpy as np
from ase.build import bulk
from ase.units import GPa
from ase.visualize import view
from mattersim.forcefield.potential import MatterSimCalculator
from mattersim.applications.phonon import PhononWorkflow

# initialize the structure of silicon
si = bulk("Si")

# attach the calculator to the atoms object
si.calc = MatterSimCalculator()

ph = PhononWorkflow(
    atoms=si,
    find_prim = False,
    work_dir = "./tmp/phonon_si_example",
    amplitude = 0.01,
    supercell_matrix = np.diag([4,4,4]),
)

has_imag, phonons = ph.run()
print(f"Has imaginary phonon: {has_imag}")
print(f"Phonon frequencies: {phonons}")


