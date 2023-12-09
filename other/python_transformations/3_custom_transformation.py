"""
TITLE: Custom Transformation (Example).
"""
"""BLOCK: Main"""
"""
This code demonstrates how to get materials from the application runtime, apply changes using python and return back to the platform
"""
# variable materials_in holds array of materials pased from the application runtime
input_materials = materials_in

def main():
  material = materials_in[0]
  poscar = material.getAsPOSCAR()

  # your code here

  
  globals()["materials_out"] = [{"poscar":poscar}]
  return globals()

main()
