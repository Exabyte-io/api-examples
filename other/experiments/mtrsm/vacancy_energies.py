import numpy as np
from ase.build import bulk
from mattersim.forcefield.potential import MatterSimCalculator
from mattersim.applications.relax import Relaxer

# 1. Build and relax the pristine supercell
si_bulk = bulk("Si", "diamond", a=5.43).repeat((3, 3, 3))
si_bulk.positions += 0.05 * np.random.randn(len(si_bulk), 3)
si_bulk.calc = MatterSimCalculator()

relaxer = Relaxer(optimizer="BFGS", constrain_symmetry=True)
relaxer.relax(si_bulk, steps=500)  # In-place relaxation
energy_bulk = si_bulk.get_potential_energy()
n_atoms_bulk = len(si_bulk)

# 2. Create a supercell with a vacancy and relax
vacancy_idx = [n_atoms_bulk // 2]  # Remove central atom for minimal defect interaction
si_vacancy = si_bulk.copy()
del si_vacancy[vacancy_idx]
si_vacancy.calc = MatterSimCalculator()
relaxer.relax(si_vacancy, steps=500)
energy_vacancy = si_vacancy.get_potential_energy()
n_atoms_vac = len(si_vacancy)

# 3. Calculate vacancy formation energy
e_vacancy = energy_vacancy - (n_atoms_vac / n_atoms_bulk * energy_bulk)
print(f"Vacancy formation energy: {e_vacancy:.4f} eV")

