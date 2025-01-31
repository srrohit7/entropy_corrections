#!/usr/bin/env python
import numpy as np
from scipy.special import i0, iv
from numpy import exp
from ase.data import atomic_numbers, atomic_masses
from ase import Atoms, Atom
from ase.io import read,write
from ase import constraints, units
import sys

#hindered rotor partition function
#Sprowl, Campbell Paper. n_eq: no of minima during a rotation, atom_fix: index of atom around which the rotation occurs, sym_num: symmetry number of the molecule 

rot_bar_eV=float(sys.argv[1])
n_eq=float(sys.argv[2])
atom_fix=int(sys.argv[3])
sym_num = float(sys.argv[4])
T= float(sys.argv[5])

g=atom_fix
slab=read("CONTCAR",format='vasp')

#The molecule rotates around fixed axis perpendicular to the surface. Hence the moment of inertia has to be calculated from the x,y coordinates of the molecule with the origin fixed at the adsorbed atom.
x_fix, y_fix, z_fix =slab.positions[g]

atom_mass=slab.get_masses()

[I_x, I_y]=[0,0]

for i in range(0,len(atom_mass)):
	if i not in slab.constraints[0].index:
		x,y,z=slab.positions[i]
		I_x+=atom_mass[i]*(x-x_fix)*(x-x_fix)
		I_y+=atom_mass[i]*(y-y_fix)*(y-y_fix)

I_red = (I_x + I_y)* units._amu * 10**-20

rot_bar= rot_bar_eV * units._e

#print(I_red)

#Rotational frequency calculated from the barrier. This differs from the value obtained in the vibrational run.
rot_freq=(0.5/np.pi)*(0.5*n_eq**2*rot_bar/I_red)**0.5

r_r = rot_bar / (units._hplanck * rot_freq)
T_r = units._k * T / (units._hplanck * rot_freq)
r_T = r_r / T_r

#print(T_r)
#qxy_num = (np.pi * r_T * exp(-r_T) * exp(-1/T_r) * (iv(0,(r_T/2))**2 )**2)**0.5
qr_num = (np.pi * r_T * exp(-r_T) * exp(-1/T_r))**0.5 * (i0(r_T/2)) * exp(1.0/((2.0+16.0*r_r)*T_r))
qr_den = (1-exp(-1/T_r))

qr= qr_num / (qr_den * sym_num) #symmetry-number?


S_rot = 8.314 * np.log(qr) #standard entropy

print(S_rot)
#print(rot_freq/(3*10**10))
#print(I_x + I_y)

#free-rotor
theta = (8 * np.pi **2 * units._c * I_red  / units._hplanck )
qr_free_1 = (np.pi * theta) ** 0.5 * (units._k * T / (units._hplanck * units._c)) ** (1.0 / 2.0)
qr_free = qr_free_1 / (sym_num * n_eq) #I'm not really sure of this 

S_rot_free = 8.314 * np.log(qr_free)
print(S_rot_free)

num_of_arg = len(sys.argv)

g = open('Hind-rot-entropy','a')

for i in range(0,num_of_arg-1):
        g.write('{}\n'.format(sys.argv[i+1]))

g.write(str(S_rot))
g.write(str('\n'))

g.write(str(S_rot_free))
g.write(str('\n'))
