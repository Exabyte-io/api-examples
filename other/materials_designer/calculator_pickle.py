import matgl
from matgl.ext.ase import M3GNetCalculator
from matgl.ext.ase import M3GNetCalculator

import cloudpickle

pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")

# Save the calculator to pickle for use in the Pyodide environment
calculator = M3GNetCalculator(pot)
print(calculator)

with open("m3gnet_calculator_3.11_32.pkl", "wb") as f:
    cloudpickle.dump(calculator, f, protocol=3)
    