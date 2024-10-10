#!/usr/bin/env python3

# vim: ts=2 expandtab ai

import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import Draw
from rdkit.Chem.rdDetermineBonds import DetermineBonds

if len(sys.argv) != 2:
  raise ValueError(f'usage: {sys.argv[0]} input')

atypes = []
coords = []
energy = False
gradient = False
xyz = False
charge = 0.

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
      charge = float(l[5:])

mol = Chem.RWMol()

idx = [ mol.AddAtom(Chem.Atom(a)) for a in atypes ]

conf = Chem.Conformer(mol.GetNumAtoms())
for i,c in enumerate(coords):
  conf.SetAtomPosition(i,Point3D(*c))

mol.AddConformer(conf)

DetermineBonds(mol)

kcalmol_in_Eh = 0.00159362
bohr_in_A = 0.529177210903

if energy:
  en = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy() * kcalmol_in_Eh
  print('energy', en)
  print()
#  print('*** energy', en,file=sys.stderr)

if gradient:
  print('gradient')
  ff = AllChem.UFFGetMoleculeForceField(mol)
  grad = np.array(ff.CalcGrad()).reshape((-1,3)) * kcalmol_in_Eh / bohr_in_A
  for i,g in enumerate(grad):
    print(i,*g)
#    print('***',i,*g,file=sys.stderr)
  print()
  
  print('charges and spin populations') # XXX: fuckup
  for i in range(len(atypes)):
    print(i,charge/len(atypes))
  print()
