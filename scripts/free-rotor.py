#!/usr/bin/env python
import numpy as np
from scipy.special import i0, iv
from numpy import exp
from ase.data import atomic_numbers, atomic_masses
from ase import Atoms, Atom
from ase.io import read,write
from ase import constraints, units
import sys

sym_num = float(sys.argv[1])
T= float(sys.argv[2])

slab=read("CONTCAR",format='vasp')

atom_mass=slab.get_masses()

m_x , m_y , m_z, M = 0., 0., 0., 0.

for i in range(0,len(atom_mass)):
        if i not in slab.constraints[0].index:
		x,y,z = slab.positions[i]
		m_x += atom_mass[i] * x
		m_y += atom_mass[i] * y 
		m_z += atom_mass[i] * z
		M += atom_mass[i]

xcm, ycm, zcm = m_x/M , m_y/M, m_z/M

Ixx, Iyy, Izz, Ixy, Ixz, Iyz = 0., 0., 0., 0., 0., 0.

for i in range(0,len(atom_mass)):
        if i not in slab.constraints[0].index:
		m = atom_mass[i]
		x_a,y_a,z_a = slab.positions[i]
		x = x_a - xcm
        	y = y_a - ycm
        	z = z_a - zcm
        	Ixx += m * (y**2. + z**2.)
        	Iyy += m * (x**2. + z**2.)
        	Izz += m * (x**2. + y**2.)
        	Ixy += m * x * y
        	Ixz += m * x * z
        	Iyz += m * y * z
# Create the inertia tensor in the current frame of reference.
I_ = np.matrix([[ Ixx, -Ixy, -Ixz],
               [-Ixy,  Iyy, -Iyz],
               [-Ixz, -Iyz,  Izz]])
# Find the eigenvalues, which are the principle moments of inertia.
I = np.linalg.eigvals(I_)

I_x , I_y , I_z = I * units._amu * 10**-20

q_x = (8 * np.pi **2 * units._c * I_x  / units._hplanck ) 
q_y = (8 * np.pi **2 * units._c * I_y  / units._hplanck ) 
q_z = (8 * np.pi **2 * units._c * I_z  / units._hplanck )

#Rotational entropy for all three rotations
q_rot_1 = (np.pi * q_x * q_y * q_z) ** 0.5 * (units._k * T / (units._hplanck * units._c)) ** (3.0 / 2.0)
q_rot = q_rot_1 / sym_num

S_rot = 8.314 * np.log (q_rot)

#Rotational entropy of each rotation
q_rot_x = ((np.pi)**(1.0/3.0) * q_x) ** 0.5 * (units._k * T / (units._hplanck * units._c)) ** (1.0 / 2.0)

S_rot_x = 8.314 * np.log (q_rot_x)

q_rot_y = ((np.pi)**(1.0/3.0) * q_y) ** 0.5 * (units._k * T / (units._hplanck * units._c)) ** (1.0 / 2.0)

S_rot_y = 8.314 * np.log (q_rot_y)

q_rot_z = ((np.pi)**(1.0/3.0) * q_z) ** 0.5 * (units._k * T / (units._hplanck * units._c)) ** (1.0 / 2.0)

S_rot_z = 8.314 * np.log (q_rot_z)

print(S_rot_x, S_rot_y, S_rot_z)
print(S_rot)

num_of_arg = len(sys.argv)

g = open('free-rot-entropy','a')
for i in range(0,num_of_arg-1):
	g.write('{}\n'.format(sys.argv[i+1]))
g.write(str(S_rot))
g.write(str('\n'))
