# -*- coding: utf-8 -*-
"""
Created on 10:06 13-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: hse06.py
"""
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import numpy as np
from matplotlib import pyplot as plt
import math
import initial
import band as bd

# ************************************ hse06 ************************************
# get k number
def hse_get_kpt_bandnum(filename = 'KPOINTS'):
    kpoints = open(filename, 'r').read().splitlines()
    num_all_kpoints = int(kpoints[1].split()[0])
    num_kpt = 0
    for flag_klines in range(3, len(kpoints)):
        if int(kpoints[flag_klines].split()[3]) == 0:
            num_kpt = flag_klines
            break
    return int(num_all_kpoints), int(num_all_kpoints-num_kpt+3)

# get the total number of high-symmetry points
def hse_get_num_high_symmetry_points(filename = "KPOINTS"):
    high_kpt_num = 0
    high_kpt_point = []
    kpt = []
    kpath = []
    kpoints = open(filename, 'r').read().splitlines()
    num_kpt_line = int(kpoints[1].split()[0])
    for flag_kpath in range(4, len(kpoints), 3):
        kpt_data1 = list(map(float, kpoints[flag_kpath].split()[0:3]))
        kpt_points_data = kpoints[flag_kpath].split('!')[1].strip()
        if kpt_points_data == "\Gamma":              # transform the gamma character
            kpt_points_data = u"Γ"
        high_kpt_num = high_kpt_num + 1              # get high symmetry point number
        high_kpt_point.append(kpt_points_data)       # get high symmetry points
        kpath.append(kpt_data1[0:3])
    kpt_points_data1 = kpoints[-1].split('!')[1].strip()
    if kpt_points_data1 == "\Gamma":             # transform the gamma character
        kpt_points_data1 = u"Γ"
    high_kpt_point.append(kpt_points_data1)      # get high symmetry points
    kpt_data2 = list(map(float, kpoints[-1].split()[0:3]))
    kpath.append(kpt_data2)                          # get high symmetry path
    high_kpt_num = int(high_kpt_num + 1)                 # get high symmetry point number
    return high_kpt_num, high_kpt_point

# get the length of high-symmetry points
def hse_get_kpt_coordinate(filename = "KPOINTS"):
    kpoints = open(filename, 'r')
    kpoints_lines = kpoints.readlines()
    num_high_symmetry, char_high_symmetry = hse_get_num_high_symmetry_points("KPATH.in")
    num_all_k, num_high_k = hse_get_kpt_bandnum()
    high_k_coordinates = np.zeros(shape=(num_high_k, 3))
    for flag_num_high in range(0, num_high_k):
        for flag_axis in range(0, 3):
            high_k_coordinates[flag_num_high][flag_axis] = kpoints_lines[num_all_k - num_high_k + 3 + flag_num_high].split()[flag_axis]
    # tranform the KPOINTS in the real coordination
    high_k_coordinate = np.zeros(shape=(num_high_k, 3))
    for i in range(0, num_high_k):
        high_k_coordinate[i][0] = high_k_coordinates[i][0] * initial.get_reciprocal()[0] + high_k_coordinates[i][1] * \
                                  initial.get_reciprocal()[3] + high_k_coordinates[i][2]*initial.get_reciprocal()[6]
        high_k_coordinate[i][1] = high_k_coordinates[i][0] * initial.get_reciprocal()[1] + high_k_coordinates[i][1] * \
                                  initial.get_reciprocal()[4] + high_k_coordinates[i][2]*initial.get_reciprocal()[7]
        high_k_coordinate[i][2] = high_k_coordinates[i][0] * initial.get_reciprocal()[2] + high_k_coordinates[i][1] * \
                                  initial.get_reciprocal()[5] + high_k_coordinates[i][2]*initial.get_reciprocal()[8]
    return high_k_coordinate

# calculate the k-length
def hse_cal_klength(high_k_coordinate):
    leng = 0
    leng_pre = 0
    length = np.zeros(shape=(len(high_k_coordinate)))
    t, g, k = 0, 0, 0
    for flag_num in range(0, len(high_k_coordinate)):
        leng = np.sqrt((high_k_coordinate[flag_num][0] - t)**2 + (high_k_coordinate[flag_num][1] - g)**2
                       + (high_k_coordinate[flag_num][2]- k)**2)
        t = high_k_coordinate[flag_num][0]
        g = high_k_coordinate[flag_num][1]
        k = high_k_coordinate[flag_num][2]
        length[flag_num] = leng + leng_pre
        leng_pre = length[flag_num]
    return length

def hse_get_eigenval():
    num_kpt_all, num_kpt_hk = hse_get_kpt_bandnum()
    num_high_k_num, high_kpt_point = hse_get_num_high_symmetry_points("KPATH.in")
    num_remove = num_kpt_all - num_kpt_hk
    fermi = initial.get_fermiLevel()
    nbd = initial.get_NBANDS() + 2
    ispin = initial.get_ispin()
    eigenval_lines = open('EIGENVAL', 'r').read().splitlines()
    if ispin == 1:
        energy_kpt = np.zeros(shape=(num_kpt_hk, nbd - 2))
        for flag_index_k_num in range(0, num_kpt_hk):
            for flag_index_orbital in range(0, nbd-2):
                energy_kpt[flag_index_k_num][flag_index_orbital] = \
                    float(eigenval_lines[int(8 + (num_remove + flag_index_k_num) * nbd + flag_index_orbital)].split()[1])-fermi
        return energy_kpt
    elif ispin == 2:
        energy_kpt_up = np.zeros(shape=(num_kpt_hk, nbd - 2))
        energy_kpt_dw = np.zeros(shape=(num_kpt_hk, nbd - 2))
        for flag_index_k_num in range(0, num_kpt_hk):
            for flag_index_orbital in range(0, nbd-2):
                energy_kpt_up[flag_index_k_num][flag_index_orbital] = \
                    float(eigenval_lines[int(8 + (num_remove + flag_index_k_num) * nbd + flag_index_orbital)].split()[1])-fermi
                energy_kpt_dw[flag_index_k_num][flag_index_orbital] = \
                    float(eigenval_lines[int(8 + (num_remove + flag_index_k_num) * nbd + flag_index_orbital)].split()[2])-fermi
        return energy_kpt_up, energy_kpt_dw
    else:
        return "No Ispin Information Obtained !"

# output band data for hse06
def hse_output_band():
    bd.output_klines("KPATH.in")
    num_kpt_all, num_kpt_hk = hse_get_kpt_bandnum()
    num_high_k_num, high_kpt_point = hse_get_num_high_symmetry_points("KPATH.in")
    kpt = hse_cal_klength(hse_get_kpt_coordinate())
    num_remove = num_kpt_all - num_kpt_hk
    fermi = initial.get_fermiLevel()
    nbd = initial.get_NBANDS()
    ispin = initial.get_ispin()
    if ispin == 1:
        energy_kpt = hse_get_eigenval()
        with open("band.dat", "w", encoding='utf-8') as band:
            band.write("Band Data for HSE" + "\n")
            for k_num in range(0, num_kpt_hk):
                band.write(str(kpt[k_num])[:8].ljust(8, '0') + "  ")
                for flag_nbd in range(0, nbd):
                    band.write(str(energy_kpt[k_num][flag_nbd])[:10].ljust(10, '0') + "  ")
                band.write("\n")
        band.close()
    elif ispin ==2:
        energy_kpt_up, energy_kpt_dw = hse_get_eigenval()
        with open("band.dat", "w", encoding='utf-8') as band:
            band.write("Band Data for HSE" + "\n")
            for k_num in range(0, num_kpt_hk):
                band.write(str(kpt[k_num])[:8].ljust(8, '0') + "  ")
                for flag_nbd in range(0, nbd):
                    band.write(str(energy_kpt_up[k_num][flag_nbd])[:10].ljust(10, '0') + "  ")
                    band.write(str(energy_kpt_dw[k_num][flag_nbd])[:10].ljust(10, '0') + "  ")
                band.write("\n")
        band.close()
        print("Wait !")
    else:
        print("Wait !")


# plot band data
def hse_bandplot(eng_r = [-10, 6, 2], color_noispin = "black", color_ispin = "red"):
    hse_output_band()
    ispin = initial.get_ispin()
    num_kpt, length, high_kpt_point, kpt = bd.get_kpt("KPATH.in")
    data_band = np.loadtxt("band.dat", skiprows=1, dtype=np.float64)
    nbd = initial.get_NBANDS()
    if ispin == 1:
        len_high_k = 0
        for j in range(0, len(length) - 1):
            length[j] = length[j] + len_high_k
            len_high_k = length[j]
        fig, ax = plt.subplots()
        for flag_nbd in range(1, nbd+1):
            plt.plot(data_band[:, 0], data_band[:, flag_nbd], color= color_noispin)
        plt.xlim(0, max(length))
        x_group_label = []
        x_start = 0
        x_group_label.append(x_start)
        for k in range(0, len(length) - 1):
            x_group_label.append(length[k])
        plt.xticks(x_group_label)
        ax.set_xticklabels(high_kpt_point, rotation=0, fontsize=12, fontname='arial')
        # plt.yticks(eng_r[0], eng_r[1], eng_r[2])
        plt.ylim(eng_r[0], eng_r[1])
        plt.ylabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=14, fontname='arial')
        ytick = np.arange(eng_r[0], eng_r[1] + 1, eng_r[2])
        a = int(len(ytick) / 2)
        plt.yticks(np.insert(ytick, a, 0), fontsize=12, )
        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=0.5, color='0.5')
        for i in length[0:-1]:
            ax.axvline(x=i, ymin=0, ymax=1, linestyle='--', linewidth=0.5, color='0.5')
        plt.savefig("band.png", dpi=300)
        plt.show()

    elif ispin == 2:
        len_high_k = 0
        for j in range(0, len(length) - 1):
            length[j] = length[j] + len_high_k
            len_high_k = length[j]
        fig, ax = plt.subplots()
        for flag_nbd in range(1, nbd * 2 + 1, 2):
            plt.plot(data_band[:, 0], data_band[:, flag_nbd], color=color_noispin)
        for flag_nbd in range(2, nbd * 2 + 1, 2):
            plt.plot(data_band[:, 0], data_band[:, flag_nbd], color=color_ispin)
        plt.xlim(0, max(length))
        x_group_label = []
        x_start = 0
        x_group_label.append(x_start)
        for k in range(0, len(length) - 1):
            x_group_label.append(length[k])
        plt.xticks(x_group_label)
        ax.set_xticklabels(high_kpt_point, rotation=0, fontsize=12, fontname='arial')
        ax.legend(("Spin_UP", "Spin_DW"), loc='best')
        # plt.yticks(eng_r[0], eng_r[1], eng_r[2])
        plt.ylim(eng_r[0], eng_r[1])
        plt.ylabel("E - " + "E$_\mathrm{{F}}$" + " (eV)", fontsize=14, fontname='arial')
        ytick = np.arange(eng_r[0], eng_r[1] + 1, eng_r[2])
        a = int(len(ytick) / 2)
        plt.yticks(np.insert(ytick, a, 0), fontsize=12, )
        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=0.5, color='0.5')
        for i in length[0:-1]:
            ax.axvline(x=i, ymin=0, ymax=1, linestyle='--', linewidth=0.5, color='0.5')
        plt.savefig("band.png", dpi=300)
        plt.show()
    else:
        print("No Ispin Information Obtained !")



def manipulate_bandplot():
    print("    ****************************************************************")
    print("    *This is a code used to plot kinds of band structure,written by*")
    print("    *                          XY Ding                             *")
    print("    ****************************************************************")
    print("\n")
    print("                       (^o^)GOOD LUCK!(^o^)                         ")
    print("\n")
    ispin = initial.get_ispin()
    if ispin == 1:
        print("    *************************** ISPIN is 1 ****************************")
        print("(1) energy range")
        print("(2) color")
        print("(3) Use default setting")
        inint = int(input("Input number:"))
        if inint == 1:
            energy = []
            min_energy = float(input("minimum energy:"))
            max_energy = float(input("maximum energy:"))
            scale_energy = float(input("energy scale:"))
            energy = [min_energy, max_energy, scale_energy]
            print("(1) continue to setting color")
            print("(2) end setting")
            con_end = int(input("Input number:"))
            if con_end == 1:
                color = str(input("which color do you want:"))
                hse_bandplot(energy, color)
            else:
                hse_bandplot(energy, "black")
        elif inint == 2:
            color = str(input("which color do you want:"))
            print("(1) continue to setting energy range")
            print("(2) end setting")
            con_end = int(input("Input number:"))
            if con_end == 1:
                min_energy = float(input("minimum energy:"))
                max_energy = float(input("maximum energy:"))
                scale_energy = float(input("energy scale:"))
                energy = [min_energy, max_energy, scale_energy]
                hse_bandplot(energy, color)
            else:
                hse_bandplot([-10, 6, 2], color)
        else:
            hse_bandplot([-10, 6, 2], "black")
    elif ispin == 2:
        print("    *************************** ISPIN is 2 ****************************")
        print("(1) energy range")
        print("(2) color")
        print("(3) Use default setting")
        inint = int(input("Input number:"))
        if inint == 1:
            energy = []
            min_energy = float(input("minimum energy:"))
            max_energy = float(input("maximum energy:"))
            scale_energy = float(input("energy scale:"))
            energy = [min_energy, max_energy, scale_energy]
            print("(1) continue to setting color")
            print("(2) end setting")
            con_end = int(input("Input number:"))
            if con_end == 1:
                color_up = str(input("color for Spin_UP:"))
                color_dw = str(input("color for Spin_DW:"))
                hse_bandplot(energy, color_up, color_dw)
            else:
                hse_bandplot(energy, "black", "red")
        elif inint == 2:
            color_up = str(input("color for Spin_UP:"))
            color_dw = str(input("color for Spin_DW:"))
            print("(1) continue to setting energy range")
            print("(2) end setting")
            con_end = int(input("Input number:"))
            if con_end == 1:
                min_energy = float(input("minimum energy:"))
                max_energy = float(input("maximum energy:"))
                scale_energy = float(input("energy scale:"))
                energy = [min_energy, max_energy, scale_energy]
                hse_bandplot(energy, color_up, color_dw)
            else:
                hse_bandplot([-10, 6, 2], color_up, color_dw)
        else:
            hse_bandplot([-10, 6, 2], "black", "red")
    else:
        print("    ******************** ISPIN is not 1 or 2!  Good bye ! ********************")

if __name__ == '__main__':
    manipulate_bandplot()
