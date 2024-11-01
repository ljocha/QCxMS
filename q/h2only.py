#!/usr/bin/env python3

# vim: ts=2 expandtab ai

import sys
import numpy as np
from copy import deepcopy


from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit

from qiskit.primitives import BackendEstimator as Estimator
from qiskit import QuantumCircuit,transpile
from qiskit_aer import AerSimulator, Aer
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import L_BFGS_B, SLSQP,COBYLA
from qiskit_algorithms import NumPyMinimumEigensolver


from qiskit_nature.second_q.mappers import BravyiKitaevMapper,BravyiKitaevSuperFastMapper,JordanWignerMapper,ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC, UCCSD


# Settings
###################################
basis = 'sto-3g'
unit = DistanceUnit.ANGSTROM
qubit_reduction = False
map_type = "BravyiKitaev"
shots = 200
optimizer = COBYLA()
ini_params = [ 0.01373469, -0.02186433, -0.30696973]
dx = 1e-05  # displacement for numerical gradient
###################################

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

        # Add displacement dx to the atom's coordinates
        self.atoms[index]["coordinates"] = [ self.atoms[index]["coordinates"][i] + dx[i] for i in range(3)]

    def displace_atom_coordinate(self, atom_index: int, axis: int, dx: float):
        if atom_index < 0 or atom_index >= len(self.atoms):
            raise IndexError("Atom index is out of range.")
        if axis not in (0, 1, 2):
            raise ValueError("Axis must be 0 (x), 1 (y), or 2 (z).")

        # Update the specified coordinate
        self.atoms[atom_index]["coordinates"][axis] += dx

    def __repr__(self):
        atom_repr = "\n".join(f"{atom['type']}: {atom['coordinates']}" for atom in self.atoms)
        return f"Energy: {self.energy}\nGradient: {self.gradient}\nAtoms:\n{atom_repr}"

def compute_energy_and_grad(mol,dx, shots, qubit_reduction, mapper_type, parameters=None, opt=SLSQP()):


    # Molecule displacement setup
    #############################################################################
    
    # r
    molecule = mol
    problem = molecule.get_pyscf_driver(unit=unit,basis=basis)
    e_nr = problem.nuclear_repulsion_energy
    second_q_ops = problem.second_q_ops()
    main_op = second_q_ops[0]

    # r - dx
    molecule_m = deepcopy(mol)
    # change - dont know the displacement r - dx
    molecule_m.displace_atom(1, (0., 0., -dx))  # XXX: H2
    #molecule_m.displace_atom_coordinate(1, 2, -dx)
    problem_mdx = molecule_m.get_pyscf_driver(unit=unit,basis=basis)
    e_nr_mdx = problem_mdx.nuclear_repulsion_energy
    second_q_ops_mdx = problem_mdx.second_q_ops()
    main_op_mdx = second_q_ops_mdx[0]

    # r + dx
    molecule_p = deepcopy(mol)
    # change - dont know the displacement r + dx
    molecule_p.displace_atom(1, (0., 0., dx)) # XXX: H2
    #molecule_p.displace_atom_coordinate(1, 2, dx)
    
    problem_pdx = molecule_p.get_pyscf_driver(unit=unit,basis=basis)
    e_nr_pdx = problem_pdx.nuclear_repulsion_energy
    second_q_ops_pdx = problem_pdx.second_q_ops()
    main_op_pdx = second_q_ops_pdx[0]
    
    num_particles = (problem.num_alpha,
                    problem.num_beta)

    num_spatial_orbitals = problem.num_spatial_orbitals
    num_spin_orbitals = 2 * problem.num_spatial_orbitals
    #############################################################################

    
    #############################################################################
    mapper = None
    if mapper_type == 'JordanWigner':
        mapper = JordanWignerMapper()
    elif mapper_type == 'ParityMapper':
        if two_qubit_reduction == False:
            mapper = ParityMapper(num_particles=problem.num_particles,two_qubit_reduction=False)
        else:
            mapper = ParityMapper(num_particles=problem.num_particles,two_qubit_reduction=True)
    elif mapper_type == 'BravyiKitaev':
        mapper = BravyiKitaevMapper()
    elif mappper_type == 'Superfast':
        mapper = BravyiKitaevSuperFastMapper(num_particles=problem.num_particles)

    # JordanWigner - https://web.archive.org/web/20191103083720/http://michaelnielsen.org/blog/archive/notes/fermions_and_jordan_wigner.pdf
    # Parity - https://arxiv.org/pdf/1701.08213
    # BravyiKitaev - https://www.sciencedirect.com/science/article/abs/pii/S0003491602962548
    # Superfast - https://arxiv.org/pdf/1712.00446

    
    # map to qubit operators
    qubit_op = mapper.map(main_op)
    qubit_op_mdx = mapper.map(main_op_mdx)
    qubit_op_pdx = mapper.map(main_op_pdx)
    #############################################################################

    # Classical - Exact solution with NUMPY
    #############################################################################
    numpy_solver = NumPyMinimumEigensolver()
    ret_exact = numpy_solver.compute_minimum_eigenvalue(qubit_op)
    #ret_exact_mdx = numpy_solver.compute_minimum_eigenvalue(qubit_op_mdx)
    #ret_exact_pdx = numpy_solver.compute_minimum_eigenvalue(qubit_op_pdx)
    e_fci = ret_exact._eigenvalue.real + e_nr
    #e_fci_mdx = ret_exact_mdx.eigenvalue.real + e_nr_mdx
    #e_fci_pdx = ret_exact_pdx.eigenvalue.real + e_nr_pdx
    print("Exact energy:", e_fci)
    #grad_fci = (e_fci_pdx - e_fci_mdx) / 2 / dx
    #############################################################################

    # PREPARE QUANTUM CIRCUIT
    #############################################################################
    

    #######
    # Simulator
    aer_sim = AerSimulator()
    # Quantum Backend
    #qc = transpile(circuit,aer_sim)
    #passmanager = generate_preset_pass_manager(optimization_level=1, backend=aer_sim)
    #qc = passmanager.run(circuit)
    estimator = Estimator(backend=aer_sim)
    ######

    # QURI functionality - maybe overwrite
    # https://docs.quantum.ibm.com/guides/qunasys-quri-chemistry#quri-chemistry-with-uccsd


    ansatz = UCC(num_spatial_orbitals=num_spatial_orbitals, num_particles=num_particles, qubit_mapper=mapper,
                excitations='sd')

    algorithm = VQE(estimator, ansatz, optimizer=opt, initial_point=parameters)

    #############################################################################

    # COMPUTE - QUANTUM
    #############################################################################
    result_cs = algorithm.compute_minimum_eigenvalue(qubit_op, aux_operators=[qubit_op_mdx, qubit_op_pdx])
    #result_mdx = algorithm.compute_minimum_eigenvalue(qubit_op_mdx)
    #result_pdx = algorithm.compute_minimum_eigenvalue(qubit_op_pdx)

    e = result_cs.eigenvalue + e_nr
    e_pdx = float(result_cs.aux_operators_evaluated[1][0]) + e_nr_pdx
    e_mdx = float(result_cs.aux_operators_evaluated[0][0]) + e_nr_mdx

    grad_cs = (e_pdx - e_mdx) / 2 / dx


    #grad_brute = (result_pdx.eigenvalue + e_nr_pdx - result_mdx.eigenvalue - e_nr_mdx) / 2 / dx
    #############################################################################

    print("Energy", e)
    print("Gradient (CS):", grad_cs)
    #print("Gradient (brute):", grad_brute)
    #print("Gradient (FCI):", grad_fci)

    # Energy, energy + dx, energy - dx, gradient correlated_sampling # parameters of VQE
    return e, e_pdx, e_mdx, grad_cs, result_cs.optimal_point



if len(sys.argv) != 2:
  raise ValueError(f'usage: {sys.argv[0]} input')

molecule = Molecule()




atypes = []
coords = []
energy = False
gradient = False
xyz = False
charge = 0

with open(sys.argv[1]) as inp:
  for l in inp:
    l=l.strip()

    if l == '*': xyz = False

    if xyz:
      s = l.split()
      atypes.append(s[0])
      coords.append(list(map(float,s[1:])))

    if l == 'energy': energy = True
    if l == 'gradient': gradient = True

    if l[:5] == '* xyz':
      xyz = True
      charge = int(l[5:])

coords = np.array(coords)
r = np.linalg.norm(coords[1,:]-coords[0,:])

# XXX: H2
molecule.add_atom('H',0.,0.,0.)
molecule.add_atom('H',0.,0.,r)


kcalmol_in_Eh = 0.00159362
bohr_in_A = 0.529177210903

en, energy_p, energy_m, g, parameters = compute_energy_and_grad(molecule ,dx,\
                                                                           shots, qubit_reduction, \
                                                                           map_type, parameters=None, \
                                                                           opt=SLSQP())

# compute grad from g by chainrule
# dE/dx = dE/dr * dr/dx
# r = \sqrt (x1-x0)^2 + ...
# dr/dx0 = 



if energy:
  en = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy() * kcalmol_in_Eh
  print('energy', en)
  print()
  print('*** energy', en,file=sys.stderr)

if gradient:
  print('gradient')
  ff = AllChem.UFFGetMoleculeForceField(mol)
  grad = np.array(ff.CalcGrad()).reshape((-1,3)) * kcalmol_in_Eh / bohr_in_A
  for i,g in enumerate(grad):
    print(i,*g)
    print('***',i,*g,file=sys.stderr)
  print()
  
  print('charges and spin populations') # XXX: fuckup
  for i in range(len(atypes)):
    print(i,charge/len(atypes))
  print()
