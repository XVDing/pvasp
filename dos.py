# -*- coding: utf-8 -*-
"""
Created on 10:06 13-11-2020 

@author: Xian-Yong Ding
mail to: dxy_vasp@163.com
python3: dos.py
"""
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import numpy as np
from matplotlib import pyplot as plt
import math
import initial

# *************************************** DOS *******************************************
# get element information from POSCAR
def _get_element(filename="POSCAR"):
    poscar = open(filename, 'r')
    pos_e = poscar.readlines()
    elements = pos_e[5].lstrip().split()
    element = []*len(elements)
    ele = []
    num_ele = []
    numbers = pos_e[6].lstrip().split()
    for flag_element in range(0, len(elements)):
        ele = elements[flag_element]
        element.append(ele)

    for flag_element_num in range(0, len(numbers)):
        ele = int(numbers[flag_element_num])
        num_ele.append(ele)
    return element, num_ele

# get data of dos from DOSCAR
def _pdos_atom():
    element, num_ele = _get_element()
    doscar = open('DOSCAR', 'r')
    dos = doscar.readlines()
    line_6 = dos[5].lstrip().split()
    eng_max = float(line_6[0])
    eng_min = float(line_6[1])
    nedos = int(line_6[2])
    fermi_eng = float(line_6[3])
    ispin = initial.get_ispin()
    if ispin == 1:
        judge_f = dos[nedos + 7]  # judge whether have the f orbitals.
        if len(judge_f.split()) > 20:
            print("Exist f orbital, This script can't manipulate the f orbital !")
        else:
            len_all_atoms = 0
            for atom_type in range(0, len(element)):
                len_all_atoms += num_ele[atom_type]
            data_dos_storage= np.zeros(shape=(len_all_atoms, nedos, 9))
            data_dos_storage_x = np.zeros(shape=(nedos))
            for flag_x in range(0, nedos):
                data_dos_storage_x[flag_x] = dos[6 + flag_x].split()[0]

            for atom_num in range(0, len_all_atoms):
                for flag_nedos in range(0, nedos):
                    for atom_orbital in range(1, 10):
                        data_dos_storage[atom_num][flag_nedos][atom_orbital-1] = \
                        dos[flag_nedos + (nedos + 1) * (atom_num + 1) + 6].split()[atom_orbital]

            data_dos_storage_orbitals = np.zeros(shape=(len(element), nedos, 9))
            num_flag_up = 0
            for flag_data_orbitals in range(0, len(element)):
                data_dos_storage_orbitals[flag_data_orbitals] = np.sum(
                    data_dos_storage[num_flag_up:num_flag_up + num_ele[flag_data_orbitals]], axis=0)
                num_flag_up += num_ele[flag_data_orbitals]


            data_dos_storage_orbitals_tot = np.zeros(shape=(len(element), nedos, 3))
            for flag_dos_tot in range(0, len(element)):
                data_dos_storage_orbitals_tot[flag_dos_tot, :, 0] = data_dos_storage_orbitals[
                                                                           flag_dos_tot, :, 0]
                data_dos_storage_orbitals_tot[flag_dos_tot, :, 1] = data_dos_storage_orbitals[
                                                                           flag_dos_tot, :,
                                                                           1] + data_dos_storage_orbitals[
                                                                                flag_dos_tot, :, 2] + \
                                                                           data_dos_storage_orbitals[
                                                                           flag_dos_tot, :, 3]
                data_dos_storage_orbitals_tot[flag_dos_tot, :, 2] = data_dos_storage_orbitals[
                                                                           flag_dos_tot, :,
                                                                           4] + data_dos_storage_orbitals[
                                                                                flag_dos_tot, :, 5] + \
                                                                           data_dos_storage_orbitals[
                                                                           flag_dos_tot, :,
                                                                           6] + data_dos_storage_orbitals[
                                                                                flag_dos_tot, :,
                                                                                7] + data_dos_storage_orbitals[
                                                                                     flag_dos_tot, :, 8]

            data_dos_storage_orbitals_tot = np.around(data_dos_storage_orbitals_tot, decimals=6)
        for ele in range(0, len(element)):
            with open("Pdos_" + str(element[ele]) + ".dat", "w", encoding='utf-8') as pdos:
                pdos.write("Pdos_" + str(element[ele]) + "\n")
                pdos.write("Energy   s   p   d" + "\n")
                for flag_line in range(0, nedos):
                    pdos.write(str(np.round(data_dos_storage_x[flag_line] - fermi_eng, 8))[0:8].ljust(8, ' ') + "   ")
                    for flag_data_orbitals in range(0, 3):
                        pdos.write(
                            str(data_dos_storage_orbitals_tot[ele][flag_line][flag_data_orbitals]).ljust(8, " ") + "   ")
                    pdos.write("\n")
            pdos.close()
    elif ispin ==2:
        judge_f = dos[nedos + 7]   # judge whether have the f orbitals.
        if len(judge_f.split()) > 20:
            print("Exist f orbital, This script can't manipulate the f orbital !")
        else:
            len_all_atoms = 0
            for atom_type in range(0, len(element)):
                len_all_atoms += num_ele[atom_type]
            data_dos_storage_spinUP=np.zeros(shape=(len_all_atoms, nedos, 9))
            data_dos_storage_spinDW=np.zeros(shape=(len_all_atoms, nedos, 9))
            data_dos_storage_x = np.zeros(shape=(nedos))
            for flag_x in range(0, nedos):
                data_dos_storage_x[flag_x] = dos[6+flag_x].split()[0]

            for atom_num in range(0, len_all_atoms):
                for flag_nedos in range(0, nedos):
                    for atom_orbital in range(1, 19, 2):
                        data_dos_storage_spinUP[atom_num][flag_nedos][int((atom_orbital-1)/2)]=dos[flag_nedos + (nedos+1)*(atom_num + 1) + 6].split()[atom_orbital]
                    for atom_orbital in range(2, 19, 2):
                        data_dos_storage_spinDW[atom_num][flag_nedos][int(atom_orbital/2-1)]=dos[flag_nedos + (nedos+1)*(atom_num + 1) + 6].split()[atom_orbital]

            data_dos_storage_orbitals_spinUP = np.zeros(shape=(len(element), nedos, 9))
            data_dos_storage_orbitals_spinDW = np.zeros(shape=(len(element), nedos, 9))
            num_flag_up = 0
            for flag_data_orbitals in range(0, len(element)):
                data_dos_storage_orbitals_spinUP[flag_data_orbitals] = np.sum(
                    data_dos_storage_spinUP[num_flag_up:num_flag_up+num_ele[flag_data_orbitals]], axis=0)
                num_flag_up += num_ele[flag_data_orbitals]
            num_flag_dw = 0
            for flag_data_orbitals in range(0, len(element)):
                data_dos_storage_orbitals_spinDW[flag_data_orbitals] = np.sum(
                    data_dos_storage_spinDW[num_flag_dw:num_flag_dw+num_ele[flag_data_orbitals]], axis=0)
                num_flag_dw += num_ele[flag_data_orbitals]

            data_dos_storage_orbitals_spinUP_tot = np.zeros(shape=(len(element), nedos, 3))
            data_dos_storage_orbitals_spinDW_tot = np.zeros(shape=(len(element), nedos, 3))
            for flag_dos_tot in range(0, len(element)):
                data_dos_storage_orbitals_spinUP_tot[flag_dos_tot, :, 0] =data_dos_storage_orbitals_spinUP[flag_dos_tot, :, 0]
                data_dos_storage_orbitals_spinUP_tot[flag_dos_tot, :, 1] = data_dos_storage_orbitals_spinUP[
                                                                           flag_dos_tot, :, 1] + data_dos_storage_orbitals_spinUP[flag_dos_tot, :, 2] + \
                                                                           data_dos_storage_orbitals_spinUP[flag_dos_tot, :, 3]
                data_dos_storage_orbitals_spinUP_tot[flag_dos_tot, :, 2] = data_dos_storage_orbitals_spinUP[
                                                                           flag_dos_tot, :, 4] + data_dos_storage_orbitals_spinUP[flag_dos_tot, :, 5] + \
                                                                           data_dos_storage_orbitals_spinUP[flag_dos_tot, :, 6] + data_dos_storage_orbitals_spinUP[flag_dos_tot, :, 7] + data_dos_storage_orbitals_spinUP[flag_dos_tot, :, 8]
                data_dos_storage_orbitals_spinDW_tot[flag_dos_tot, :, 0] =data_dos_storage_orbitals_spinDW[flag_dos_tot, :, 0]
                data_dos_storage_orbitals_spinDW_tot[flag_dos_tot, :, 1] = data_dos_storage_orbitals_spinDW[
                                                                           flag_dos_tot, :, 1] + data_dos_storage_orbitals_spinDW[flag_dos_tot, :, 2] + \
                                                                           data_dos_storage_orbitals_spinDW[flag_dos_tot, :, 3]
                data_dos_storage_orbitals_spinDW_tot[flag_dos_tot, :, 2] = data_dos_storage_orbitals_spinDW[
                                                                           flag_dos_tot, :, 4] + data_dos_storage_orbitals_spinDW[flag_dos_tot, :, 5] + \
                                                                           data_dos_storage_orbitals_spinDW[flag_dos_tot, :, 6] + data_dos_storage_orbitals_spinDW[flag_dos_tot, :, 7] + data_dos_storage_orbitals_spinDW[flag_dos_tot, :, 8]


            data_dos_storage_orbitals_spinUP_tot = np.around(data_dos_storage_orbitals_spinUP_tot, decimals=6)
            data_dos_storage_orbitals_spinDW_tot = np.around(data_dos_storage_orbitals_spinDW_tot, decimals=6)
        for ele in range(0, len(element)):
            with open("Pdos_" + str(element[ele]) + ".dat", "w", encoding='utf-8') as pdos:
                pdos.write("Pdos_" + str(element[ele]) +  "\n")
                pdos.write("Energy  s_up  s_dw  p_up  p_dw  d_up  d_dw" + "\n")
                for flag_line in range(0, nedos):
                    pdos.write(str(np.round(data_dos_storage_x[flag_line]-fermi_eng, 8))[0:8].ljust(8, ' ') + "   ")
                    for flag_data_orbitals in range(0, 3):
                        pdos.write(str(data_dos_storage_orbitals_spinUP_tot[ele][flag_line][flag_data_orbitals]).ljust(8, " ")+"   ")
                        pdos.write(str(-data_dos_storage_orbitals_spinDW_tot[ele][flag_line][flag_data_orbitals]).ljust(8, " ") + "   ")
                    pdos.write("\n")
            pdos.close()
    else:
        print("Can't read the ISPIN information, check whether have INCAR file or not !")

def _pdos_atom_orbital():
    element, num_ele = _get_element()
    doscar = open('DOSCAR', 'r')
    dos = doscar.readlines()
    line_6 = dos[5].lstrip().split()
    eng_max = float(line_6[0])
    eng_min = float(line_6[1])
    nedos = int(line_6[2])
    fermi_eng = float(line_6[3])
    ispin = initial.get_ispin()
    if ispin == 1:
        judge_f = dos[nedos + 7]   # judge whether have the f orbitals.
        if len(judge_f.split()) > 10:
            print("Exist f orbital, This script can't manipulate the f orbital !")
        else:
            len_all_atoms = 0
            for atom_type in range(0, len(element)):
                len_all_atoms += num_ele[atom_type]
            data_dos_storage=np.zeros(shape=(len_all_atoms, nedos, 9))
            data_dos_storage_x = np.zeros(shape=(nedos))
            for flag_x in range(0, nedos):
                data_dos_storage_x[flag_x] = dos[6+flag_x].split()[0]

            for atom_num in range(0, len_all_atoms):
                for flag_nedos in range(0, nedos):
                    for atom_orbital in range(1, 10):
                        data_dos_storage[atom_num][flag_nedos][atom_orbital-1]=dos[flag_nedos + (nedos+1)*(atom_num + 1) + 6].split()[atom_orbital]
            data_dos_storage_orbitals = np.zeros(shape=(len(element), nedos, 9))
            num_flag = 0
            for flag_data_orbitals in range(0, len(element)):
                data_dos_storage_orbitals[flag_data_orbitals] = np.sum(
                    data_dos_storage[num_flag:num_flag+num_ele[flag_data_orbitals]], axis=0)
                num_flag += num_ele[flag_data_orbitals]
            np.set_printoptions(precision=4)
            data_dos_storage_orbitals = np.around(data_dos_storage_orbitals, decimals=6)
        for ele in range(0, len(element)):
            with open("Pdos_" + str(element[ele]) + ".dat", "w", encoding='utf-8') as pdos:
                pdos.write("Pdos_" + str(element[ele]) +  "\n")
                pdos.write("Energy  s  py  pz  px  dxy  dyz  dz2  dxz  dx2" + "\n")
                for flag_line in range(0, nedos):
                    pdos.write(str(np.round(data_dos_storage_x[flag_line]-fermi_eng, 8)).ljust(8, ' ') + "   ")
                    for flag_data_orbitals in range(0, 9):
                        pdos.write(str(data_dos_storage_orbitals[ele][flag_line][flag_data_orbitals]).ljust(8, " ")+"   ")
                    pdos.write("\n")
            pdos.close()

    elif ispin ==2:
        judge_f = dos[nedos + 7]   # judge whether have the f orbitals.
        if len(judge_f.split()) > 20:
            print("Exist f orbital, This script can't manipulate the f orbital !")
        else:
            len_all_atoms = 0
            for atom_type in range(0, len(element)):
                len_all_atoms += num_ele[atom_type]
            data_dos_storage_spinUP=np.zeros(shape=(len_all_atoms, nedos, 9))
            data_dos_storage_spinDW=np.zeros(shape=(len_all_atoms, nedos, 9))
            data_dos_storage_x = np.zeros(shape=(nedos))
            for flag_x in range(0, nedos):
                data_dos_storage_x[flag_x] = dos[6+flag_x].split()[0]
            print(len(dos[6010].split()))
            for atom_num in range(0, len_all_atoms):
                for flag_nedos in range(0, nedos):
                    for atom_orbital in range(1, 19, 2):
                        data_dos_storage_spinUP[atom_num][flag_nedos][int((atom_orbital-1)/2)]=dos[flag_nedos + (nedos+1)*(atom_num + 1) + 6].split()[atom_orbital]
                    for atom_orbital in range(2, 19, 2):
                        data_dos_storage_spinDW[atom_num][flag_nedos][int(atom_orbital/2-1)]=dos[flag_nedos + (nedos+1)*(atom_num + 1) + 6].split()[atom_orbital]

            data_dos_storage_orbitals_spinUP = np.zeros(shape=(len(element), nedos, 9))
            data_dos_storage_orbitals_spinDW = np.zeros(shape=(len(element), nedos, 9))
            num_flag = 0
            for flag_data_orbitals in range(0, len(element)):
                data_dos_storage_orbitals_spinUP[flag_data_orbitals] = np.sum(
                    data_dos_storage_spinUP[num_flag:num_flag+num_ele[flag_data_orbitals]], axis=0)
                num_flag += num_ele[flag_data_orbitals]
            num_flag = 0
            for flag_data_orbitals in range(0, len(element)):
                data_dos_storage_orbitals_spinDW[flag_data_orbitals] = np.sum(
                    data_dos_storage_spinDW[num_flag:num_flag+num_ele[flag_data_orbitals]], axis=0)
                num_flag += num_ele[flag_data_orbitals]

            np.set_printoptions(precision=4)
            data_dos_storage_orbitals_spinUP = np.around(data_dos_storage_orbitals_spinUP, decimals=6)
            data_dos_storage_orbitals_spinDW = np.around(data_dos_storage_orbitals_spinDW, decimals=6)
        for ele in range(0, len(element)):
            with open("Pdos_" + str(element[ele]) + ".dat", "w", encoding='utf-8') as pdos:
                pdos.write("Pdos_" + str(element[ele]) +  "\n")
                pdos.write("Energy  s_up  s_dw  py_up  py_dw  pz_up  pz_dw  px_up  px_dw  dxy_up  dxy_dw  dyz_up  dyz_dw  dz2_up  dz2_dw  dxz_up  dxz_dw  dx2_up  dx2_dw" + "\n")
                for flag_line in range(0, nedos):
                    pdos.write(str(np.round(data_dos_storage_x[flag_line]-fermi_eng, 8))[0:8].ljust(8, ' ') + "   ")
                    for flag_data_orbitals in range(0, 9):
                        pdos.write(str(data_dos_storage_orbitals_spinUP[ele][flag_line][flag_data_orbitals]).ljust(8, " ")+"   ")
                        pdos.write(str(-data_dos_storage_orbitals_spinDW[ele][flag_line][flag_data_orbitals]).ljust(8, " ") + "   ")
                    pdos.write("\n")
            pdos.close()

# get data of dos from DOSCAR
def _dos_tot():
    element, num_ele = _get_element()
    doscar = open('DOSCAR', 'r')
    dos = doscar.readlines()
    line_6 = dos[5].lstrip().split()
    eng_max = float(line_6[0])
    eng_min = float(line_6[1])
    nedos = int(line_6[2])
    fermi_eng = float(line_6[3])
    print(fermi_eng)
    ispin = initial.get_ispin()
    if ispin == 1:
        with open("Tdos.dat", "w", encoding='utf-8') as tdos:
            tdos.write("Total density of state" + "\n")
            tdos.write("Energy".ljust(13," ") + "TDOS".ljust(13," ") + "IDOS".ljust(13," ") + "\n")
            for flag_line in range(0, nedos):
                tdos.write(str(float(dos[6+flag_line].split()[0])-fermi_eng)[:10].ljust(10, " ") + "   " +
                           str(dos[6+flag_line].split()[1]) + "   "+
                           str(dos[6+flag_line].split()[2]) + "\n")
        tdos.close()
    elif ispin == 2:
        with open("Tdos_SpinUp.dat", "w", encoding='utf-8') as tdos_up:
            tdos_up.write("Tdos_spinUP" + "\n")
            tdos_up.write("Energy".ljust(13," ") + "TDOS".ljust(13," ") + "IDOS".ljust(13," ") + "\n")
            for flag_line in range(0, nedos):
                tdos_up.write(str(float(dos[6+flag_line].split()[0])-fermi_eng)[:10].ljust(10, " ") + "   " +
                           str(dos[6+flag_line].split()[1])+ "   "+
                           str(dos[6+flag_line].split()[3])+ "\n")
        tdos_up.close()
        with open("Tdos_SpinDw.dat", "w", encoding='utf-8') as tdos_dw:
            tdos_dw.write("Tdos_spinDW" + "\n")
            tdos_dw.write("Energy".ljust(13," ") + "TDOS".ljust(13," ") + "IDOS".ljust(13," ") + "\n")
            for flag_line in range(0, nedos):
                tdos_dw.write(str(float(dos[6+flag_line].split()[0])-fermi_eng)[:10].ljust(10, " ") + "   " +
                           str("-"+dos[6+flag_line].split()[2]) + "   "+
                           str("-"+dos[6+flag_line].split()[4]) + "\n")
        tdos_dw.close()
    else:
        print("No ispin information obtained !")



def _dos_plot_atom_orbital(eng_range = [-10, 10, 2]):
    _pdos_atom_orbital()
    orbitals_noispin = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"]
    orbitals_ispin = ["s_up", "s_dw", "py_up", "py_dw", "pz_up", "pz_dw", "px_up", "px_dw", "dxy_up", "dxy_dw", "dyz_up", "dyz_dw", "dz2_up", "dz2_dw", "dxz_up", "dxz_dw", "dx2_up", "dx2_dw"]
    ispin = initial.get_ispin()
    element, num_ele = _get_element()
    filename = [[] for i in range(0, len(element))]
    data = [[] for i in range(0, len(element))]
    for flag_atom_type in range(0, len(element)):
        filename[flag_atom_type] = "Pdos_" + element[flag_atom_type] + ".dat"
        data[flag_atom_type] = np.loadtxt(filename[flag_atom_type], skiprows=2, dtype=float)

    get_d_maxvalue = [[] for i in range(0, len(element))]
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    if ispin == 1:
        print(" *************************** ISPIN is 1 ****************************")
        for flag_atom_type in range(0, len(element)):
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 5]+data[flag_atom_type][:, 6]+data[flag_atom_type][:, 7]+data[flag_atom_type][:, 8]+data[flag_atom_type][:, 9])
            if get_d_maxvalue[flag_atom_type] == 0:     # if d orbitals have no values, plot the s p
                for flag_atom_orbital in range(1, 5):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital], ls="-",
                             label=str(element[flag_atom_type] +"_"+ orbitals_noispin[flag_atom_orbital-1]))
            else:
                for flag_atom_orbital in range(1, 10):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital], ls="-",
                             label=str(element[flag_atom_type] +"_"+ orbitals_noispin[flag_atom_orbital-1]))

        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
        ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
        plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
        plt.xlim(eng_range[0], eng_range[1])
        xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
        plt.xticks(xtick)

        plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
        plt.yticks([])
        plt.legend(fontsize=12, loc="best")
        plt.savefig("Pdos_atom_orbitals.png", dpi=300)
        plt.show()
    elif ispin == 2:
        print(" *************************** ISPIN is 2 ****************************")
        for flag_atom_type in range(0, len(element)):
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 10]+data[flag_atom_type][:, 11]+data[flag_atom_type][:, 12]+data[flag_atom_type][:, 13]+data[flag_atom_type][:, 14])
            if get_d_maxvalue[flag_atom_type] == 0:     # if d orbitals have no values, plot the s p
                for flag_atom_orbital in range(1, 9):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital], ls="-",
                             label=str(element[flag_atom_type] +"_"+ orbitals_ispin[flag_atom_orbital-1]))
            else:
                for flag_atom_orbital in range(1, 19):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital], ls="-",
                             label=str(element[flag_atom_type] + "_" + orbitals_ispin[flag_atom_orbital-1]))

        ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
        plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
        plt.xlim(eng_range[0], eng_range[1])
        xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
        plt.xticks(xtick)

        plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
        plt.yticks([])
        plt.legend(fontsize=12, loc="best")
        plt.savefig("Pdos_atom_orbitals.png", dpi=300)
        plt.show()

    else:
        print("No Ispin Information Obtained !, Checking whether this path exist the INCAR file")


def _dos_plot_atom(eng_range = [-10, 10, 2]):
    _pdos_atom()
    orbitals_noispin = ["s", "p", "d"]
    orbitals_ispin = ["s_up", "s_dw", "p_up", "p_dw", "d_up", "d_dw"]
    ispin = initial.get_ispin()
    element, num_ele = _get_element()
    filename = [[] for i in range(0, len(element))]
    data = [[] for i in range(0, len(element))]
    for flag_atom_type in range(0, len(element)):
        filename[flag_atom_type] = "Pdos_" + element[flag_atom_type] + ".dat"
        data[flag_atom_type] = np.loadtxt(filename[flag_atom_type], skiprows=2, dtype=float)

    get_d_maxvalue = [[] for i in range(0, len(element))]
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    if ispin == 1:
        print("    *************************** ISPIN is 1 ****************************")
        for flag_atom_type in range(0, len(element)):
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 3])

            if get_d_maxvalue[flag_atom_type] == 0:  # if d orbitals have no values, plot the s p
                for flag_atom_orbital in range(1, 3):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital],
                             ls="-",label=str(element[flag_atom_type] + "_" + orbitals_noispin[flag_atom_orbital-1]))
            else:
                for flag_atom_orbital in range(1, 4):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital],
                             ls="-",label=str(element[flag_atom_type] + "_" + orbitals_noispin[flag_atom_orbital-1]))


        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
        ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
        plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
        plt.xlim(eng_range[0], eng_range[1])
        xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
        plt.xticks(xtick)

        plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
        plt.yticks([])
        plt.legend(fontsize=12, loc="best")
        plt.savefig("Pdos_atom.png", dpi=300)
        plt.show()
    elif ispin == 2:
        print("    *************************** ISPIN is 2 ****************************")
        for flag_atom_type in range(0, len(element)):
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 5]+data[flag_atom_type][:, 6])

            if get_d_maxvalue[flag_atom_type] == 0:  # if d orbitals have no values, plot the s p
                for flag_atom_orbital in range(1, 5):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital],
                             ls="-",label=str(element[flag_atom_type] + "_" + orbitals_ispin[flag_atom_orbital-1]))
            else:
                for flag_atom_orbital in range(1, 7):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, flag_atom_orbital],
                             ls="-",label=str(element[flag_atom_type] + "_" + orbitals_ispin[flag_atom_orbital-1]))


        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
        ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
        plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
        plt.xlim(eng_range[0], eng_range[1])
        xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
        plt.xticks(xtick)

        plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
        plt.yticks([])
        plt.legend(fontsize=12, loc="best")
        plt.savefig("Pdos_atom.png", dpi=300)
        plt.show()

    else:
        print("No Ispin Information Obtained !, Checking whether this path exist the INCAR file")

def _dos_plot_tot(eng_range = [-10, 10, 2], integral = 0):
    _dos_tot()
    ispin = initial.get_ispin()
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    if integral == 0:
        if ispin == 1:
            data = np.loadtxt("Tdos.dat", skiprows=2, dtype=float)
            plt.plot(data[:, 0], data[:, 1], ls="-",
                     label="TDOS")
            ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
            ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
            plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
            plt.xlim(eng_range[0], eng_range[1])
            xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
            plt.xticks(xtick)

            plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
            plt.yticks([])
            plt.legend(fontsize=12, loc="best")
            plt.savefig("Total_Tdos.png", dpi=300)
            plt.show()

        elif ispin == 2:
            data_up = np.loadtxt("Tdos_SpinUp.dat", skiprows=2, dtype=float)
            data_dw = np.loadtxt("Tdos_SpinDw.dat", skiprows=2, dtype=float)
            plt.plot(data_up[:, 0], data_up[:, 1], ls="-",
                     label="TDOS_SpinUP")
            plt.plot(data_dw[:, 0], data_dw[:, 1], ls="-",
                     label="TDOS_SpinDW")
            ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
            ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
            plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
            plt.xlim(eng_range[0], eng_range[1])
            xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
            plt.xticks(xtick)

            plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
            plt.yticks([])
            plt.legend(fontsize=12, loc="best")
            plt.savefig("Total_Tdos.png", dpi=300)
            plt.show()

        else:
            print("No Ispin Information Obtained !, Checking whether this path exist the INCAR file")
    else:
        if ispin == 1:
            data = np.loadtxt("Tdos.dat", skiprows=2, dtype=float)
            plt.plot(data[:, 0], data[:, 1], ls="-",
                     label="TDOS")
            plt.plot(data[:, 0], data[:, 2], ls="-",
                     label="IDOS")
            ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
            ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
            plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
            plt.xlim(eng_range[0], eng_range[1])
            xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
            plt.xticks(xtick)

            plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
            plt.yticks([])
            plt.legend(fontsize=12, loc="best")
            plt.savefig("Total_TIdos.png", dpi=300)
            plt.show()

        elif ispin == 2:
            data_up = np.loadtxt("Tdos_SpinUp.dat", skiprows=2, dtype=float)
            data_dw = np.loadtxt("Tdos_SpinDw.dat", skiprows=2, dtype=float)
            plt.plot(data_up[:, 0], data_up[:, 1], ls="-",
                     label="TDOS_SpinUP")
            plt.plot(data_up[:, 0], data_up[:, 2], ls="-",
                     label="IDOS_SpinUP")

            plt.plot(data_dw[:, 0], data_dw[:, 1], ls="-",
                     label="TDOS_SpinDW")
            plt.plot(data_dw[:, 0], data_dw[:, 2], ls="-",
                     label="IDOS_SpinDW")
            ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
            ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
            plt.xlabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=16, fontname='arial')
            plt.xlim(eng_range[0], eng_range[1])
            xtick = np.arange(eng_range[0], eng_range[1] + 1, eng_range[2])
            plt.xticks(xtick)

            plt.ylabel("DOS (states/eV)", fontsize=16, fontname='arial')
            plt.yticks([])
            plt.legend(fontsize=12, loc="best")
            plt.savefig("Total_TIdos.png", dpi=300)
            plt.show()
        else:
            print("No Ispin Information Obtained !, Checking whether this path exist the INCAR file")


def _manipulate_dos():
    print(" ************************** Enter dos ploting *****************************")
    print("(1) Total dos")
    print("(2) prjected dos for different kind atoms of each orbitals (eg. s, p, d)")
    print("(3) prjected dos for different kind atoms of each partical orbitals (eg. s, py, pz ...)")
    select_dos = int(input("input a number: "))
    if select_dos == 1:
        print("(1) energy range")
        print("(2) Use default setting")
        select_eng = int(input("input a number: "))
        if select_eng == 1:
            min_eng = int(input("minimum energy is: "))
            max_eng = int(input("maximum energy is: "))
            eng_scale = float(input("energy scale is: "))
            _dos_plot_tot([min_eng, max_eng, eng_scale])
        elif select_eng == 2:
            _dos_plot_tot()
        else:
            print("please input a right number !")

    elif select_dos == 2:
        print("(1) energy range")
        print("(2) Use default setting")
        select_eng = int(input("input a number: "))
        if select_eng == 1:
            min_eng = int(input("minimum energy is: "))
            max_eng = int(input("maximum energy is: "))
            eng_scale = float(input("energy scale is: "))
            _dos_plot_atom([min_eng, max_eng, eng_scale])
        elif select_eng == 2:
            _dos_plot_atom()
        else:
            print("please input a right number !")
    elif select_dos == 3:
        print("(1) energy range")
        print("(2) Use default setting")
        select_eng = int(input("input a number: "))
        if select_eng == 1:
            min_eng = int(input("minimum energy is: "))
            max_eng = int(input("maximum energy is: "))
            eng_scale = float(input("energy scale is: "))
            _dos_plot_atom_orbital([min_eng, max_eng, eng_scale])
        elif select_eng == 2:
            _dos_plot_atom_orbital()
        else:
            print("please input a right number !")

    else:
        print("please input a right number !")

if __name__ == '__main__':
    _manipulate_dos()