#!/usr/bin/env python3

# vim: ts=2 expandtab ai

import sys
import copy
import numpy as np

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.mappers import BravyiKitaevMapper,BravyiKitaevSuperFastMapper,JordanWignerMapper,ParityMapper

from qiskit_algorithms import NumPyMinimumEigensolver


unit = DistanceUnit.ANGSTROM
basis = 'sto-3g'
dx = 1e-5

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

coords = np.array(coords)

# XXX: spin not used, qcxms gives weird numbers
def op_from_geom(atoms,xyz,charge,spin):
  geom_str = '\n'.join(
    [ atoms[i] + ' ' + ' '.join(map(str,xyz[i])) for i in range(len(atoms)) ]
  )
  
  driver = PySCFDriver(atom=geom_str, unit=unit, basis=basis, charge=charge, spin=charge) # XXX
  problem = driver.run()
  e_nr = problem.nuclear_repulsion_energy
  
  op = problem.hamiltonian.second_q_op()
  mapper = ParityMapper(num_particles=problem.num_particles)
  qubit_op = mapper.map(op)
  return qubit_op,e_nr

main_op,e_nr = op_from_geom(atypes,coords,charge,spin)

aux_ops = []

for a in range(len(atypes)):
  for c in range(3):
    d = np.zeros_like(coords)
    d[a,c] = dx
    aux_ops += [
      op_from_geom(atypes,coords-d,charge,spin)[0],
      op_from_geom(atypes,coords+d,charge,spin)[0]
      ]


numpy_solver = NumPyMinimumEigensolver()
assert numpy_solver.supports_aux_operators()

ret_exact = numpy_solver.compute_minimum_eigenvalue(main_op,aux_operators=aux_ops)

aux_ops_m = np.array([ op[0].real for op in ret_exact.aux_operators_evaluated[::2] ]).reshape(-1,3)
aux_ops_p = np.array([ op[0].real for op in ret_exact.aux_operators_evaluated[1::2] ]).reshape(-1,3)

grad = (aux_ops_p - aux_ops_m) / (2. * dx)

en = ret_exact._eigenvalue.real + e_nr

kcalmol_in_Eh = 0.00159362
bohr_in_A = 0.529177210903

if energy:
  print('energy', en * kcalmol_in_Eh)
  print()
  print('*** energy', e_nr, en, en * kcalmol_in_Eh,file=sys.stderr)

if gradient:
  print('gradient')
  grad *= kcalmol_in_Eh / bohr_in_A
  for i,g in enumerate(grad):
    print(i,*g)
    print('***',i,*g,file=sys.stderr)
  print()
  
  print('charges and spin populations') # XXX: fuckup
  for i in range(len(atypes)):
    print(i,charge/len(atypes))
  print()
