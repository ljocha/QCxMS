#!/usr/bin/env python3

# vim: ts=2 expandtab ai

import sys
import numpy as np

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.mappers import BravyiKitaevMapper,BravyiKitaevSuperFastMapper,JordanWignerMapper,ParityMapper

from qiskit_algorithms import NumPyMinimumEigensolver

unit = DistanceUnit.ANGSTROM
basis = 'sto-3g'

if len(sys.argv) != 2:
  raise ValueError(f'usage: {sys.argv[0]} input')

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
      charge,spin = [int(i) for i in l[5:].split()]

geom_str = '\n'.join(
  [ atypes[i] + ' ' + ' '.join(map(str,coords[i])) for i in range(len(atypes)) ]
)

driver = PySCFDriver(atom=geom_str, unit=unit, basis=basis, charge=charge, spin=spin)
problem = driver.run()

e_nr = problem.nuclear_repulsion_energy
main_op = problem.hamiltonian.second_q_op()
mapper = ParityMapper(num_particles=problem.num_particles)
qubit_op = mapper.map(main_op)

numpy_solver = NumPyMinimumEigensolver()
ret_exact = numpy_solver.compute_minimum_eigenvalue(qubit_op)

en = ret_exact._eigenvalue.real + e_nr

kcalmol_in_Eh = 0.00159362
bohr_in_A = 0.529177210903

if energy:
  print('energy', en * kcalmol_in_Eh)
  print()
  print('*** energy', en * kcalmol_in_Eh,file=sys.stderr)

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
