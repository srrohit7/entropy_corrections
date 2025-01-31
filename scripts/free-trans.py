#!/usr/bin/env python
import numpy as np
from scipy.special import i0
from numpy import exp 
from ase.data import atomic_numbers, atomic_masses
from ase import Atoms, Atom
from ase.io import read,write
from ase import constraints, units
import sys
from ase import units

nn_dis_1=float(sys.argv[1])
nn_dis_2=float(sys.argv[2])
sym=sys.argv[3]
T= float(sys.argv[4])

slab=read("CONTCAR",format='vasp')

atom_mass=slab.get_masses()

#print(slab.get_masses())
ad_mass=0

for i in range(0,len(atom_mass)):
	if i not in slab.constraints[0].index:
		ad_mass += atom_mass[i]			

red_mass = ad_mass * units._amu
b_1=nn_dis_1 * 10**-10
b_2=nn_dis_2 * 10**-10

if sym == 'tri':
        q_ideal= (0.5 * 3**0.5 * b_1 * b_2) * (2 * np.pi * red_mass * units._k * T) / (units._hplanck**2)
if sym == 'rect':
        q_ideal= (b_1 * b_2) * (2 * np.pi * red_mass * units._k * T) / (units._hplanck**2)


#	q_ideal= (0.5 * 3**0.5 * b**2) * (2 * np.pi * red_mass * units._k * T) / (units._hplanck**2)

S_ideal = 8.314 * np.log(q_ideal)


S_ideal_meV = S_ideal / 96

print(S_ideal_meV*298.16/1000)

num_of_arg = len(sys.argv)

g = open('Free-trans-entropy','a')
for i in range(0,num_of_arg-1):
        g.write('{}\n'.format(sys.argv[i+1]))
g.write(str(S_ideal_meV))
g.write(str('\n'))
