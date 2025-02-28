### Qiskit Nature
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import TwoLocal


from qiskit_nature.second_q.mappers import BravyiKitaevMapper,BravyiKitaevSuperFastMapper,JordanWignerMapper,ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC, UCCSD
from qiskit_nature.second_q.algorithms import GroundStateEigensolver, ExcitedStatesEigensolver, QEOM, EvaluationRule
from qiskit_nature.second_q.transformers import FreezeCoreTransformer


### Qiskit Runtime
from qiskit.primitives import BackendEstimator as Estimator
from qiskit import QuantumCircuit,transpile
from qiskit_aer import AerSimulator, Aer
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService


### Qiskit Algorithms
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import L_BFGS_B, SLSQP,COBYLA,SPSA
from qiskit_algorithms import NumPyMinimumEigensolver


from qiskit_nature.second_q.algorithms.initial_points import HFInitialPoint

### Other
import numpy as np
import copy
import warnings
warnings.filterwarnings("ignore")
from decimal import Decimal

from molecule import Molecule


# Settings
###################################
basis = 'sto-3g'
unit = DistanceUnit.ANGSTROM
qubit_reduction = False
map_type = "JordanWigner"
shots = 1000
#optimizer = SLSQP(maxiter=10,ftol=1e-9))
optimizer = SPSA(maxiter=10)
ini_params = [ 0, 0, 0]
dx = 1e-5  # displacement for numerical gradient
###################################


def compute_energy_and_grad(mol,dx, shots, qubit_reduction, mapper_type, parameters=None, opt=SLSQP()):


    # Molecule displacement setup
    #############################################################################
    
    # r
    molecule = mol
    problem = molecule.get_pyscf_driver(unit=unit,basis=basis)
    e_nr = problem.nuclear_repulsion_energy
    main_op = problem.hamiltonian.second_q_op()
    
    molecule.print_atoms()
    
    # r - dx
    molecule_m = copy.deepcopy(mol)
    molecule_m.displace_atom(1, (0, 0, -dx))
    molecule_m.print_atoms()
    #molecule_m.displace_atom_coordinate(1, 2, -dx)
    problem_mdx = molecule_m.get_pyscf_driver(unit=unit,basis=basis)
    e_nr_mdx = problem_mdx.nuclear_repulsion_energy
    second_q_ops_mdx = problem_mdx.second_q_ops()
    main_op_mdx = problem_mdx.hamiltonian.second_q_op()

    # r + dx
    molecule_p = copy.deepcopy(mol)
    molecule_p.displace_atom(1, (0, 0, dx))
    molecule_p.print_atoms()
    #molecule_p.displace_atom_coordinate(1, 2, dx)
    problem_pdx = molecule_p.get_pyscf_driver(unit=unit,basis=basis)
    e_nr_pdx = problem_pdx.nuclear_repulsion_energy
    second_q_ops_pdx = problem_pdx.second_q_ops()
    main_op_pdx = problem_pdx.hamiltonian.second_q_op()
    
    
    transformer = FreezeCoreTransformer()
    transformed_problem = transformer.transform(problem)
    transformed_problem_mdx = transformer.transform(problem_mdx)
    transformed_problem_pdx = transformer.transform(problem_pdx)
    
    # GET OPERATORS
    #############################################################################

    num_particles = (problem.num_alpha,
                    problem.num_beta)

    num_spatial_orbitals = problem.num_spatial_orbitals
    num_spin_orbitals = 2 * problem.num_spatial_orbitals
    print("\n\n")
    print(f'Nuclear repulsion energy: {problem.nuclear_repulsion_energy} Ha')
    print(f'Reference energy: {problem.reference_energy} Ha')
    print(f'Number of spin orbitals: {problem.num_spin_orbitals}')
    print(f'Number of alpha electrons: {problem.num_alpha}')
    print(f'Number of beta electrons: {problem.num_beta}')
    print("\n\n")
    #############################################################################
    
    #############################################################################
    
    print(f'Nuclear repulsion energy: {transformed_problem.nuclear_repulsion_energy} Ha')
    print(f'Reference energy: {transformed_problem.reference_energy} Ha')
    print(f'Number of spin orbitals: {transformed_problem.num_spin_orbitals}')
    print(f'Number of alpha electrons: {transformed_problem.num_alpha}')
    print(f'Number of beta electrons: {transformed_problem.num_beta}')
    #############################################################################
    
    # MAPPERS
    #############################################################################
    mapper = None
    if mapper_type == 'JordanWigner':
        mapper = JordanWignerMapper()
    elif mapper_type == 'ParityMapper':
        if qubit_reduction == False:
            mapper = ParityMapper(num_particles=transformed_problem.num_particles)
        else:
            mapper = ParityMapper(num_particles=transformed_problem.num_particles,two_qubit_reduction=True)
    elif mapper_type == 'BravyiKitaev':
        mapper = BravyiKitaevMapper()
    elif mappper_type == 'Superfast':
        mapper = BravyiKitaevSuperFastMapper(num_particles=transformed_problem.num_particles)

    # JordanWigner - https://web.archive.org/web/20191103083720/http://michaelnielsen.org/blog/archive/notes/fermions_and_jordan_wigner.pdf
    # Parity - https://arxiv.org/pdf/1701.08213
    # BravyiKitaev - https://www.sciencedirect.com/science/article/abs/pii/S0003491602962548
    # Superfast - https://arxiv.org/pdf/1712.00446

    # map to qubit operators
    tapered_mapper = transformed_problem.get_tapered_mapper(mapper)
    tapered_mapper_mdx = transformed_problem_mdx.get_tapered_mapper(mapper)
    tapered_mapper_pdx = transformed_problem_pdx.get_tapered_mapper(mapper)
    
    qubit_op = tapered_mapper.map(transformed_problem.hamiltonian.second_q_op())
    qubit_op_mdx = tapered_mapper_mdx.map(transformed_problem_mdx.hamiltonian.second_q_op())
    qubit_op_pdx = tapered_mapper_pdx.map(transformed_problem_pdx.hamiltonian.second_q_op())
    
    #print(qubit_op)
    #print(qubit_op_mdx)
    #print(qubit_op_pdx)
    
    # PREPARE QUANTUM CIRCUIT
    #############################################################################
    
    # BACKEND
    #############################################################################
    
    aer_sim = AerSimulator(method="statevector")
    estimator = Estimator(backend=aer_sim)
    
    ansatz = UCCSD(
        transformed_problem.num_spatial_orbitals,
        transformed_problem.num_particles,
        tapered_mapper,
        initial_state=HartreeFock(
            transformed_problem.num_spatial_orbitals,
            transformed_problem.num_particles,
            tapered_mapper,
        ),
    )

    vqe_solver = VQE(estimator, ansatz=ansatz, optimizer=opt)
    

    
    # SETUP INITIAL POINT
    #############################################################################
    
    vqe_solver.initial_point = np.zeros(ansatz.num_parameters)
    
    initial_point = HFInitialPoint()
    initial_point.ansatz = ansatz
    initial_point.problem = transformed_problem
    vqe_solver.initial_point = initial_point.to_numpy_array()
    
    #############################################################################
    
    solver = GroundStateEigensolver(tapered_mapper, vqe_solver)
    op, _ = solver.get_qubit_operators(transformed_problem)
    print("\n\n")
    print(f'Number of qubits of main op: {op.num_qubits}, number of paulis: {len(op.paulis)}')
    
    
    print(tapered_mapper_mdx.pauli_table)
    solver1 = GroundStateEigensolver(tapered_mapper_mdx, vqe_solver)
    print(solver1)
    op1, _ = solver.get_qubit_operators(transformed_problem_mdx)
    print("\n\n")
    print(op1)
    print(f'Number of qubits of - op: {op1.num_qubits}, number of paulis: {len(op1.paulis)}')
    
    solver2 = GroundStateEigensolver(tapered_mapper_pdx, vqe_solver)
    op2, _ = solver.get_qubit_operators(transformed_problem_pdx)
    print("\n\n")
    print(f'Number of qubits of + op: {op2.num_qubits}, number of paulis: {len(op2.paulis)}')
    
    
    result_cs = []
    if solver.supports_aux_operators():
        result_cs = solver.solve(transformed_problem, aux_operators={'qubit_op_mdx': op1, 'qubit_op_pdx': op2})
    print(result_cs)
    #result_cs = solver.solve(transformed_problem)
    #result_cs2 = solver1.solve(transformed_problem_mdx)
    #result_cs3 = solver2.solve(transformed_problem_pdx)
    
    # Classical - Exact solution with NUMPY
    #############################################################################
    numpy_solver = NumPyMinimumEigensolver()
    ret_exact = numpy_solver.compute_minimum_eigenvalue(op)
    ret_exact_mdx = numpy_solver.compute_minimum_eigenvalue(op1)
    ret_exact_pdx = numpy_solver.compute_minimum_eigenvalue(op2)
    e_fci = ret_exact._eigenvalue.real + e_nr
    e_fci_mdx = ret_exact_mdx.eigenvalue.real + e_nr_mdx
    e_fci_pdx = ret_exact_pdx.eigenvalue.real + e_nr_pdx
    grad_fci = (e_fci_pdx - e_fci_mdx) / 2 / dx
    #############################################################################
    
    ###############################################################################
    e = result_cs.eigenvalues + e_nr
    
    #e_pdx = result_cs2.eigenvalues + e_nr_pdx
    #e_mdx = result_cs3.eigenvalues + e_nr_mdx
    
    e_pdx = float(result_cs.aux_operators_evaluated[0].get('qubit_op_pdx')) + e_nr_pdx
    e_mdx = float(result_cs.aux_operators_evaluated[0].get('qubit_op_mdx')) + e_nr_mdx
    
    e_pdx = Decimal(str(e_pdx))
    e_mdx = Decimal(str(e_mdx))
    dx = Decimal(str(dx))
    
    grad_cs = ((e_pdx - e_mdx) / (2*dx))
    

    ##############################################################################
    print("\n\n")
    print("Nuclear repulsion energy:", e_nr)
    print("Electronic groundstate energy:", result_cs.groundenergy)
    print("Overall Energy: ", e)
    print("Energy in +:", e_pdx)
    print("Energy in -:", e_mdx)
    print("Gradient (CS):", grad_cs)
    print("\n\n")
    print("Classical Simulation")
    print("Exact energy:", e_fci)
    print("Exact energy in +:", e_fci_pdx)
    print("Exact energy in -:", e_fci_mdx)
    print("Gradient (FCI):", grad_fci)
    #
    
    # Energy, energy + dx, energy - dx, gradient correlated_sampling # parameters of VQE
    return e_fci, grad_fci



def main():
    molecule = Molecule()
    molecule.load_from_file("example2.txt")
    print(molecule)

    energy, gradient = compute_energy_and_grad(molecule ,dx,shots, qubit_reduction, map_type, parameters=None, opt=optimizer)

    return energy, gradient





if __name__ == "__main__":
    main()