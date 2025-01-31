#!/usr/bin/env python
from __future__ import division
import numpy as np
from ase.data import atomic_numbers, atomic_masses
from ase import units
from numpy import exp
import sys

T=float(sys.argv[2])
cutoff_freq = float(sys.argv[1])

#OPEN and reading a file
freq=open("freq.txt")
freq_lines=freq.readlines()
freq.close()

#freq_num=freq_lines.split()
#freq_num=map(float,freq_num)
#print(freq_num)


vib_energy=np.zeros(len(freq_lines))
Ev=0
ZPE=0

hc=units._hplanck * units._c
beta= 1 / (units._k * T)
#print(units._k * T/(1.6 * 10 ** -19))
R = units._k * units._Nav

#Reference Atkins 10th edition
for i in range(0, len(freq_lines)):
	freq_num=freq_lines[i].split()
	#freq_num=map(float,freq_num)
	if float(freq_num[0])>cutoff_freq :
		vib_energy[i]=hc*100*float(freq_num[0]) #freq in cm-1 to 1/s
		beta_e = vib_energy[i] * beta #e/kT
		Ev_SI =  vib_energy[i] / (np.exp(beta_e)-1)
		Ev += Ev_SI / (1.6 * 10 ** -19) #in eV
		ZPE += vib_energy[i] / (3.2 * 10 ** -19)
#for energy in vib_energy:
	
print(Ev)
print(ZPE)

g = open('ZPE','a')

g.write('{}\n'.format(T))
g.write(str(ZPE))
g.write(str('\n'))

