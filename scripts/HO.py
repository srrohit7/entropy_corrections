#!/usr/bin/env python
from __future__ import division
import numpy as np
from ase.data import atomic_numbers, atomic_masses
from ase import units
from numpy import exp
import sys
from ase.io import read,write

T=float(sys.argv[2])
cutoff_freq = float(sys.argv[1])
num_of_arg = len(sys.argv)
#OPEN and reading a file
freq=open("freq.txt")
freq_lines=freq.readlines()
freq.close()

#freq_num=freq_lines.split()
#freq_num=map(float,freq_num)
#print(freq_num)


vib_energy=np.zeros(len(freq_lines))
Sv = 0
Sv_t = 0

hc=units._hplanck * units._c
beta= 1 / (units._k * T)
R = units._k * units._Nav

#Reference Atkins 10th edition
for i in range(0, len(freq_lines)):
    freq_num=freq_lines[i].split()
    #print(freq_num[0])
    #freq_num=map(float,freq_num)
    if float(freq_num[0])>cutoff_freq :
        vib_energy[i]=hc*100*float(freq_num[0])
        beta_e = vib_energy[i] * beta
        qv = 1 / (1-exp(-beta_e))
        Sv_t += R * ( (beta_e)/(exp(beta_e)-1) + np.log(qv) ) 
        Sv += R * ( np.log(qv) )
#for energy in vib_energy:
	
Sv_meV = Sv / 96

Sv_t_meV = Sv_t/96

print(Sv_meV*298.16/1000)

print(Sv_meV)
print(Sv_t_meV)
g = open('HO-entropy','a')
for i in range(0,num_of_arg-2):
        g.write('{}\n'.format(sys.argv[i+1]))
g.write(str(Sv_meV)) 
g.write(str('\n'))
