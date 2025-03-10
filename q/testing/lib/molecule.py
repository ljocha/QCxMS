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
        self.free_indices = []  # Indexy volných souřadnic
    
    def add_atom(self, atom_type: str, x: float, y: float, z: float):
        atom = {
            "type": atom_type,
            "coordinates": (x, y, z)
        }
        self.atoms.append(atom)
    
    def load_from_file(self, filename: str):
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        self.atoms = []
        self.free_indices = []
        
        for line in lines:
            parts = line.split()
            if len(parts) == 4:
                atom_type, x, y, z = parts[0], parts[1], parts[2], parts[3]
                coordinates = []
                for i, val in enumerate([x, y, z]):
                    if "{}" in val:
                        self.free_indices.append((len(self.atoms), i))
                        coordinates.append(0.0)  # Dočasná hodnota, bude nahrazena
                    elif "-{}" in val:
                        self.free_indices.append((len(self.atoms), i))
                        coordinates.append(0.0)  # Uložíme mínusovou souřadnici správně
                    else:
                        coordinates.append(float(val))
                self.add_atom(atom_type, *coordinates)
    
    def apply_configuration(self, values):
        if len(values) != len(self.free_indices):
            raise ValueError("Počet hodnot neodpovídá počtu volných souřadnic.")
        
        for (atom_idx, coord_idx), new_value in zip(self.free_indices, values):
            original_value = self.atoms[atom_idx]["coordinates"][coord_idx]
            self.atoms[atom_idx]["coordinates"] = (
                self.atoms[atom_idx]["coordinates"][0] if coord_idx != 0 else new_value if original_value >= 0 else -new_value,
                self.atoms[atom_idx]["coordinates"][1] if coord_idx != 1 else new_value if original_value >= 0 else -new_value,
                self.atoms[atom_idx]["coordinates"][2] if coord_idx != 2 else new_value if original_value >= 0 else -new_value,
            )
    
    def displace_atom_coordinate(self, atom_index: int, axis: int, dx: float):
        if atom_index < 0 or atom_index >= len(self.atoms):
            raise IndexError("Atom index is out of range.")
        if axis not in (0, 1, 2):
            raise ValueError("Axis must be 0 (x), 1 (y), or 2 (z).")

        self.atoms[atom_index]["coordinates"] = (
            self.atoms[atom_index]["coordinates"][0] if axis != 0 else self.atoms[atom_index]["coordinates"][0] + dx,
            self.atoms[atom_index]["coordinates"][1] if axis != 1 else self.atoms[atom_index]["coordinates"][1] + dx,
            self.atoms[atom_index]["coordinates"][2] if axis != 2 else self.atoms[atom_index]["coordinates"][2] + dx,
        )
    
    def get_pyscf_driver(self, unit: DistanceUnit = DistanceUnit.ANGSTROM, basis: str = "sto-3g"):
        """
        Vytvoří PySCF driver pro výpočet energie a gradientu molekuly.
        """
        geometry_str = "\n".join(
            f"{atom['type']} {atom['coordinates'][0]} {atom['coordinates'][1]} {atom['coordinates'][2]}"
            for atom in self.atoms
        )
        driver = PySCFDriver(atom=geometry_str, unit=unit, basis=basis)
        return driver.run()
    
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