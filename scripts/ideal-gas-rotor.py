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

I=slab.get_moments_of_inertia()

I_x , I_y , I_z = I * units._amu * 10**-20

print(I)

q_x = (8 * np.pi **2 * units._c * I_x  / units._hplanck )
q_y = (8 * np.pi **2 * units._c * I_y  / units._hplanck )
q_z = (8 * np.pi **2 * units._c * I_z  / units._hplanck )

q_rot_1 = ((np.pi * q_x * q_y * q_z) ** 0.5) * (units._k * T / (units._hplanck * units._c)) ** (1.5)
q_rot = q_rot_1 / sym_num

S_rot = 8.314 * np.log (q_rot)
S_rot_T=8.314 * (1.5+ np.log (q_rot))

S_rot_meV = S_rot/(96)
S_rot_T_meV= S_rot_T/(96)

print(S_rot_meV)
print(S_rot_T_meV)

g = open('IG-rot-entropy','a')

num_of_arg = len(sys.argv)

for i in range(0,num_of_arg-2):
        g.write('{}\n'.format(sys.argv[i+1]))

g.write(str(S_rot_meV))
g.write(str('\n'))
