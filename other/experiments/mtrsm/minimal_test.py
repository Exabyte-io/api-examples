import torch
from loguru import logger
from ase.build import bulk
from ase.units import GPa
from mattersim.forcefield import MatterSimCalculator

device = "cuda" if torch.cuda.is_available() else "cpu"
logger.info(f"Running MatterSim on {device}")

si = bulk("Si", "diamond", a=5.43)
si.calc = MatterSimCalculator(device=device)
logger.info(f"Energy (eV)                 = {si.get_potential_energy()}")
logger.info(f"Energy per atom (eV/atom)   = {si.get_potential_energy()/len(si)}")
logger.info(f"Forces of first atom (eV/A) = {si.get_forces()[0]}")
logger.info(f"Stress[0][0] (eV/A^3)       = {si.get_stress(voigt=False)[0][0]}")
logger.info(f"Stress[0][0] (GPa)          = {si.get_stress(voigt=False)[0][0] / GPa}")
