import os
import numpy as np
from itertools import product
from lib.molecule import Molecule
from decimal import Decimal, ROUND_HALF_UP

def generate_molecule_grid(template_file, config_file, steps, dx):
    """
    Generuje grid molekul posunutých v rámci definovaných stupňů volnosti.
    :param template_file: Cesta k template souboru s `{}`.
    :param config_file: Cesta ke konfiguračnímu souboru s kompletní molekulou.
    :param steps: Počet kroků posunu.
    :param dx: Velikost posunu.
    :return: Seznam objektů Molecule s posunutými geometriemi.
    """
    # Načtení molekuly přímo z config souboru
    molecule = Molecule()
    molecule.load_from_file(config_file)
    
    # Načtení template pro určení stupňů volnosti
    template_molecule = Molecule()
    template_molecule.load_from_file(template_file)
    free_indices = template_molecule.free_indices  # Indexy souřadnic, které budeme posouvat
    
    if not free_indices:
        raise ValueError("V template souboru nebyly nalezeny žádné stupně volnosti.")
    
    # Vytvoření gridu posunutí s vysokou přesností Decimal a zaokrouhlením
    displacement_values = [Decimal(str(i)) * Decimal(str(dx)) for i in range(steps)]
    displacement_combinations = list(product(displacement_values, repeat=len(free_indices)))
    
    molecules = []
    for displacements in displacement_combinations:
        displaced_molecule = Molecule()
        displaced_molecule.atoms = [atom.copy() for atom in molecule.atoms]
        displaced_molecule.free_indices = free_indices.copy()
        
        # Posouváme všechny volné souřadnice o odpovídající hodnoty displacementu
        for (atom_idx, coord_idx), disp in zip(free_indices, displacements):
            original_value = Decimal(str(molecule.atoms[atom_idx]["coordinates"][coord_idx]))
            adjusted_dx = -disp if original_value < 0 else disp
            new_value = (original_value + adjusted_dx).quantize(Decimal('1.0000'), rounding=ROUND_HALF_UP)
            displaced_molecule.displace_atom_coordinate(atom_idx, coord_idx, float(new_value - original_value))
        
        molecules.append(displaced_molecule)
    
    return molecules