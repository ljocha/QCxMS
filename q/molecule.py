from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit

import numpy as np
import copy
import warnings
warnings.filterwarnings("ignore")
from decimal import Decimal


class Molecule:
    def __init__(self):
        self.atoms = []
        self.energy = None
        self.gradient = None

    def add_atom(self, atom_type: str, x: float, y: float, z: float):
        atom = {
            "type": atom_type,
            "coordinates": (x, y, z)
        }
        self.atoms.append(atom)

    def load_from_file(self, filename: str):
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Reading the atoms section
        for line in lines[3:]:  # Skip the first three lines
            if line.strip() == "*":  # End
                break
            parts = line.split()
            atom_type = parts[0]
            x, y, z = map(float, parts[1:4])
            self.add_atom(atom_type, x, y, z)

    def get_pyscf_driver(self, unit: DistanceUnit = DistanceUnit.ANGSTROM, basis: str = "sto-3g"):
        # Create molecular geometry string
        geometry_str = "\n".join(f"{atom['type']} {atom['coordinates'][0]} {atom['coordinates'][1]} {atom['coordinates'][2]}"
                                 for atom in self.atoms)
        
        # Generate and return the PySCFDriver
        driver = PySCFDriver(atom=geometry_str, unit=unit, basis=basis)
        
        problem = driver.run()
        
        return problem
    
    def displace_atom(self, index: int, dx: tuple):
        if index < 0 or index >= len(self.atoms):
            raise IndexError("Atom index is out of range.")
    
        # Convert to Decimal for high-precision arithmetic
        self.atoms[index]["coordinates"] = [
            float(Decimal(str(self.atoms[index]["coordinates"][i])) + Decimal(str(dx[i]))) for i in range(3)]
    
    
    def displace_atom_coordinate(self, atom_index: int, axis: int, dx: float):
        if atom_index < 0 or atom_index >= len(self.atoms):
            raise IndexError("Atom index is out of range.")
        if axis not in (0, 1, 2):
            raise ValueError("Axis must be 0 (x), 1 (y), or 2 (z).")

        # Update the specified coordinate - change
        self.atoms[atom_index]["coordinates"][axis] += dx
    
    def print_atoms(self):
        if not self.atoms:
            print("No atoms in the molecule.")
        else:
            print("Atoms in the molecule:")
            for i, atom in enumerate(self.atoms):
                print(f"{i + 1}. {atom['type']} - Coordinates: {atom['coordinates']}")
        
    def __repr__(self):
        atom_repr = "\n".join(f"{atom['type']}: {atom['coordinates']}" for atom in self.atoms)
        return f"Energy: {self.energy}\nGradient: {self.gradient}\nAtoms:\n{atom_repr}"
