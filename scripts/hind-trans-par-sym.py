#!/usr/bin/env python
import numpy as np
from scipy.special import i0
from numpy import exp 
from ase.data import atomic_numbers, atomic_masses
from ase import Atoms, Atom
from ase.io import read,write
from ase import constraints, units
import sys

trans_bar_eV=float(sys.argv[1])
nn_dis=float(sys.argv[2])
T= float(sys.argv[4])
sym=sys.argv[3]
slab=read("CONTCAR",format='vasp')

atom_mass=slab.get_masses()

#print(slab.get_masses())
ad_mass=0

for i in range(0,len(atom_mass)):
	if i not in slab.constraints[0].index:
		ad_mass += atom_mass[i]			

red_mass = ad_mass * units._amu
b=nn_dis * 10**-10
trans_bar= trans_bar_eV * units._e

trans_freq=(trans_bar/(2* red_mass * b**2))**0.5

r_x = trans_bar / (units._hplanck * trans_freq)
T_x = units._k * T / (units._hplanck * trans_freq)
r_T = r_x / T_x

qxy_num = np.pi * r_T * exp(-r_T) * exp(-1/T_x) * (i0(r_T/2))**2
qxy_den = (1-exp(-1/T_x))**2

qxy= qxy_num / qxy_den

S_trans = 8.314 * np.log(qxy)
b_1 = nn_dis * 10**-10
b_2 = nn_dis * 10**-10

#q_ideal= (0.5 * 3**0.5 * b**2) * (2 * np.pi * red_mass * units._k * T) / (units._hplanck**2)
if sym == 'tri':
        q_ideal= (0.5 * 3**0.5 * b_1 * b_2) * (2 * np.pi * red_mass * units._k * T) / (units._hplanck**2)
if sym == 'rect':
        q_ideal= (b_1 * b_2) * (2 * np.pi * red_mass * units._k * T) / (units._hplanck**2)

S_ideal = 8.314 * np.log(q_ideal)
print(S_trans)
#print(trans_freq/(3*10**10))
print(S_ideal)

num_of_arg = len(sys.argv)

g = open('Hind-trans-entropy','a')
for i in range(0,num_of_arg-1):
        g.write('{}\n'.format(sys.argv[i+1]))
g.write(str(S_trans))
g.write(str('\n'))
g.write(str(S_ideal))
g.write(str('\n'))
