# -*- coding: utf-8 -*-
"""
Created on 10:14 13-11-2020 

@author: Xian-Yong Ding
mail to: dxy_vasp@163.com
python3: function.py
"""
import numpy as np
from matplotlib import pyplot as plt
import math 

# get the spin information from INCAR file
def get_ispin():
    incar = open('INCAR','r')
    incar_lines = incar.readlines()
    flag=[]
    for flag_incar in range(0,len(incar_lines)):
        if 'ISPIN' in incar_lines[flag_incar]:
            flag = incar_lines[flag_incar].lstrip()
            spinline = incar_lines[flag_incar].split('=')[1]
            if flag[0]=='#':          # whether the ISPIN is commented or not
                    ispin = 1
            else:
                    ispin = int(spinline.split()[0])
            break
        else:
                ispin = 1
    return int(ispin)


# get NBANDS from OUTCAR
def get_NBANDS():
    outcar = open('OUTCAR','r')
    outcar_lines = outcar.readlines()
    nbands = 0
    for flag_outcar in range(0,len(outcar_lines)):
        if 'NBANDS' in outcar_lines[flag_outcar]:
            if len(outcar_lines[flag_outcar].split('=')) != 4:
                continue
            else:
                spinline = outcar_lines[flag_outcar].split('=')[3]
                nbands = int(spinline.split()[0])
                break
    return nbands

# get the fermi level from OUTCAR
def get_fermiLevel():
    outcar = open('OUTCAR','r')
    outcar_lines = outcar.readlines()
    for flag_outcar in range(0,len(outcar_lines)):
        if 'E-fermi' in outcar_lines[flag_outcar]:
            spinline = outcar_lines[flag_outcar].split(':')[1]
            efermi = float(spinline.split()[0])
            break
    return efermi

# get reciprocal unit cell
def ChaCheng(a1, a2):
    c = []
    c1 = float(a1[1] * a2[2] - a1[2] * a2[1])
    c.append(c1)
    c2 = float(a1[2] * a2[0] - a1[0] * a2[2])
    c.append(c2)
    c3 = float(a1[0] * a2[1] - a1[1] * a2[0])
    c.append(c3)
    return c
def DianCheng(b1, b2):
    d1 = float(b1[0] * b2[0])
    d2 = float(b1[1] * b2[1])
    d3 = float(b1[2] * b2[2])
    d = d1 + d2 + d3
    return d

def get_poscar():
    poscar = open('POSCAR', 'r')
    poscar_lines = poscar.readlines()
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    return pos

def get_reciprocal():
    pos = get_poscar()
    rep = []
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [pos[6], pos[7], pos[8]]
    volume = DianCheng(a1, ChaCheng(a2, a3))
    scalar = 2 * math.pi /volume
    for i in range(0, 3):
        b1 = scalar*ChaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*ChaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * ChaCheng(a1, a2)[k]
        rep.append(b3)
    return rep

if __name__ == '__main__':
    print(get_poscar())
    print(get_reciprocal())