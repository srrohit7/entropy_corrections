#!/usr/bin/env python
 
import numpy as np
from ase.io import read, write
 
x = open("OUTCAR","r")
 
read = x.readlines()
x.close()
 
f = open("frequency.txt","w+")
a = []
 
for i, line in enumerate(read):
    if "2PiTHz" in line:
        for l in read[i:i+1]:
            f.write(l,)
