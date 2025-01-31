#!/usr/bin/env python
import numpy as np
from scipy.special import i0
from numpy import exp 
from ase.data import atomic_numbers, atomic_masses
from ase import Atoms, Atom
from ase.io import read,write
from ase import constraints, units
import sys

#hind-trans-par.py dif_barrier_1 dif_barrier_2 nn_dis_1 nn_dis_2 symmetry(tri or rect) T

#hindered-translator asymmetric in two directions. nn_dis is the diffusion length or the nearest neighbor distance
#Need to fix the area in the 2-D ideal gas

trans_bar_eV_1=float(sys.argv[1])
trans_bar_eV_2=float(sys.argv[2])
nn_dis_1=float(sys.argv[3])
nn_dis_2=float(sys.argv[4])
sym=sys.argv[5]
T= float(sys.argv[6])

slab=read("CONTCAR",format='vasp')

atom_mass=slab.get_masses()

#print(slab.get_masses())
ad_mass=0

for i in range(0,len(atom_mass)):
	if i not in slab.constraints[0].index:
		ad_mass += atom_mass[i]			

#print(ad_mass)

red_mass = ad_mass * units._amu
b_1=nn_dis_1 * 10**-10
b_2=nn_dis_2 * 10**-10

b= [b_1 , b_2]

trans_bar_eV = [trans_bar_eV_1 , trans_bar_eV_2]

qxy=1

for i in range(0,len(trans_bar_eV)):
	trans_bar= trans_bar_eV[i] * units._e

	trans_freq=(trans_bar/(2* red_mass * b[i]**2))**0.5
	#print(trans_freq/(3*10**10))

	r_x = trans_bar / (units._hplanck * trans_freq)
	T_x = units._k * T / (units._hplanck * trans_freq)
	r_T = r_x / T_x

	qxy_num = np.pi * r_T * exp(-r_T) * exp(-1/T_x) * (i0(r_T/2))**2 * exp(2.0/((2.0+16.0*r_x)*T_x))
	qxy_den = (1-exp(-1/T_x))**2

	qxy *= (qxy_num / qxy_den)**0.5

S_trans = 8.314 * np.log(qxy)

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
