"""TITLE: Custom Transformation (Example)."""
"""BLOCK: Main"""
"""
This code demonstrates how to:
 - get materials from the application runtime as JS class instances, 
 - extract poscar data from the materials, 
 - return the data back to the application runtime.
See https://github.com/Exabyte-io/made.js/blob/dev/src/material.js also.
"""
import json

# materials_in holds the array of materials pased from the application runtime
input_materials = materials_in


def main():
    material = input_materials[0]
    poscar = material.getAsPOSCAR()

    # REMOVE THIS COMMENT AND PLACE YOUR CODE HERE.
    # NOTE: when passing back to the application runtime, the poscar data will be converted to a JS class instance
    # of a Material class (https://github.com/Exabyte-io/made.js/blob/dev/src/material.js).

    output_materials = [{"poscar": poscar}]
    print("Output materials: ", json.dumps(output_materials, indent=4))

    # Return the globals() dictionary to the application environment. Required for proper operation at this time.
    globals()["materials_out"] = output_materials
    return globals()


main()
