import numpy as np
from ase.build import bulk
from ase.units import GPa
from mattersim.forcefield.potential import MatterSimCalculator
from mattersim.applications.relax import Relaxer

# initialize the structure of silicon
si = bulk("Si", "diamond", a=5.43)

# perturb the structure
si.positions += 0.1 * np.random.randn(len(si), 3)

# attach the calculator to the atoms object
si.calc = MatterSimCalculator()

# initialize the relaxation object
relaxer = Relaxer(
    optimizer="BFGS", # the optimization method
    filter="ExpCellFilter", # filter to apply to the cell
    constrain_symmetry=True, # whether to constrain the symmetry
)

relaxed_structure = relaxer.relax(si, steps=500)
