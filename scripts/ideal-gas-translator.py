#!/usr/bin/env python
import numpy as np
from scipy.special import i0, iv
from numpy import exp
from ase.data import atomic_numbers, atomic_masses
from ase import Atoms, Atom
from ase.io import read,write
from ase import constraints, units
import sys

T= float(sys.argv[1])

slab=read("CONTCAR",format='vasp')

atom_mass=slab.get_masses()
tot_mass=0

for i in range(0,len(atom_mass)):
	tot_mass += atom_mass[i]

mass = tot_mass * units._amu

Po = 1.01325 * 10**5

q_ideal= (units._k * T / Po) * ((2 * np.pi * mass * units._k * T) / (units._hplanck**2))**(1.5)

S_ideal = 8.314 * np.log(q_ideal)
S_ideal_T= 8.314 * (2.50 + np.log(q_ideal))

S_ideal_meV = S_ideal / 96
S_ideal_T_meV = S_ideal_T/96

g = open('ideal-gas-trans-entropy','a')
g.write('{}\n'.format(T))
g.write(str(S_ideal_meV))
g.write(str('\n'))
print(S_ideal_meV)
print(S_ideal_T_meV)

