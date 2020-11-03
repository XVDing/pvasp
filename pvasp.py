# -*- coding: utf-8 -*-
"""
Created on 9:39 30-10-2020 

@author: Xian-Yong Ding
mail to: dxy_vasp@163.com
python3: dos.py
"""
import numpy as np
from matplotlib import pyplot as plt
import math 

#------------------ FONT_setup ----------------------
font = {'family': 'arial',
        'color': 'black',
        'weight': 'normal',
        'size': 13.0,
        }
#------------------- Data manipulating ----------------------



# ******************** band ************************
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
    for flag_outcar in range(0,len(outcar_lines)):
        if 'NBANDS' in outcar_lines[flag_outcar]:
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

# get K-PATH form KPOINTS and calculating the kpt
def get_kpt():
    high_kpt_num = 0
    high_kpt_point = []
    kpt = []
    kpath = []
    kpoints = open('KPOINTS', 'r').read().splitlines()
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

    # calculating all the kpoints for the k axis of band
    # 1. transform the high symmetry points to its real coordinates in k space
    k_real = []
    k_middle = []
    k = [0, 0, 0]
    kx = [0, 0, 0]
    ky = [0, 0, 0]
    kz = [0, 0, 0]
    for flag_k in range(0, len(kpath)):
        for flag_d in range(0, 3):
            kx[flag_d] = kpath[flag_k][0] * get_reciprocal()[flag_d]
            ky[flag_d] = kpath[flag_k][1] * get_reciprocal()[flag_d+3]
            kz[flag_d] = kpath[flag_k][2] * get_reciprocal()[flag_d+6]
            k[flag_d] = kx[flag_d] + ky[flag_d] + kz[flag_d]
        k_middle = [k[0], k[1], k[2]]                   # real high symmetry point coordinations
        k_real.append(k_middle)

    # calculating the kpath_length
    length = [0] * len(k_real)
    for cal_k in range(0, len(k_real)):
        if cal_k < len(k_real)-1:
            for i in range(0, 3):
                length[cal_k] = (
                    (k_real[cal_k + 1][0] - k_real[cal_k][0]) ** 2 + (k_real[cal_k + 1][1] - k_real[cal_k][1]) ** 2 + (
                                k_real[cal_k + 1][2] - k_real[cal_k][2]) ** 2) ** (0.5)
        else:
            break

    # get all the kpoints
    num_kpt = (num_kpt_line) * (high_kpt_num-1)
    kpt = []
    kn = 0.0
    lent = 0
    for flag_kpt_hn in range(0, high_kpt_num-1):
        for flag_kpt_pn in range(0, num_kpt_line):
            kn = lent + (length[flag_kpt_hn] / (num_kpt_line -1)) * flag_kpt_pn
            round(float(kn), 6)
            kpt.append(kn)
        lent = lent + length[flag_kpt_hn]
    return int(num_kpt), length, high_kpt_point, kpt

# get Energy level form Eigenval file
def get_eigenval():
    num_kpt, length, high_kpt_point, kpt = get_kpt()
    fermi = get_fermiLevel()
    eig = []
    eig_up = []
    eig_dw = []
    energy_kpt = []
    energy_kpt_up = []
    energy_kpt_dw = []
#    energy = []
    nbd = get_NBANDS() + 2
    ispin = get_ispin()
    eigenval_lines = open('EIGENVAL', 'r').read().splitlines()
    if ispin == 1:
        for flag_nbd_num in range(0, nbd-2):
            for flag_k_num in range(0, int(num_kpt)+1):
                if flag_k_num < int(num_kpt):
                    flag_orbital = 8 + flag_nbd_num + flag_k_num * nbd
                    eig = float(eigenval_lines[flag_orbital].split()[1]) - fermi
                    round(float(eig), 6)
                    energy_kpt.append(eig)
        return energy_kpt
    elif ispin == 2:
        for flag_nbd_num in range(0, nbd-2):
            for flag_k_num in range(0, int(num_kpt)+1):
                if flag_k_num < int(num_kpt):
                    flag_orbital = 8 + flag_nbd_num + flag_k_num * nbd
                    eig_up = float(eigenval_lines[flag_orbital].split()[1]) - fermi
                    eig_dw = float(eigenval_lines[flag_orbital].split()[2]) - fermi
                    round(float(eig_up), 6)
                    round(float(eig_dw), 6)
                    energy_kpt_up.append(eig_up)
                    energy_kpt_dw.append(eig_dw)
        return energy_kpt_up, energy_kpt_dw
    else:
        return "No Ispin Information Obtained !"


# output band.dat
def output_band():
    ispin = get_ispin()
    if ispin == 1:
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = get_NBANDS()
        kpt = np.transpose(np.array(kpoints))
        # kpt = kpoints.astype(np.float)
        eng = np.array(get_eigenval())
        energy = eng.reshape(nbd, num_kpt)
        # reverse the sequence
        #    rev_kpt = kpt[::-1]
        with open("band.dat", "w", encoding='utf-8') as band:
            band.write("Band Data" + "\n")
            for k_num in range(0, num_kpt):
                band.write(str(kpt[k_num])[:10].ljust(10, '0') + "  ")
                for flag_nbd in range(0, nbd):
                    band.write(str(energy[flag_nbd][k_num])[:10].ljust(10, '0') + "  ")
                band.write("\n")
        band.close()
    elif ispin == 2:
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = get_NBANDS()
        kpt = np.transpose(np.array(kpoints))
        eng_up, eng_dw = get_eigenval()
        eig_up = np.array(eng_up)
        eig_dw = np.array(eng_dw)
        energy_up = eig_up.reshape(nbd, num_kpt)
        energy_dw = eig_dw.reshape(nbd, num_kpt)
        with open("band_up.dat", "w", encoding='utf-8') as band_up:
            band_up.write("Band Data" + "\n")
            for k_num in range(0, num_kpt):
                band_up.write(str(kpt[k_num])[:10].ljust(10, '0') + "    ")
                for flag_nbd in range(0, nbd):
                    band_up.write(str(energy_up[flag_nbd][k_num])[:10].ljust(10, '0') + "    ")
                band_up.write("\n")
        band_up.close()
        with open("band_dw.dat", "w", encoding='utf-8') as band_dw:
            band_dw.write("Band Data" + "\n")
            for k_num in range(0, num_kpt):
                band_dw.write(str(kpt[k_num])[:10].ljust(10, '0') + "    ")
                for flag_nbd in range(0, nbd):
                    band_dw.write(str(energy_dw[flag_nbd][k_num])[:10].ljust(10, '0') + "    ")
                band_dw.write("\n")
        band_dw.close()
    else:
        print("No Ispin Information Obtained !")


# plot band
def bandplot(eng_r = [-10, 6, 2], color_noispin = "black", color_ispin = "red"):
    output_band()
    ispin = get_ispin()
    if ispin == 1:
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = get_NBANDS()
        banddata = np.loadtxt("band.dat", skiprows=1, dtype=np.float64)
        kpt = banddata[:, 0]
        len_high_k = 0
        for j in range(0, len(length) - 1):
            length[j] = length[j] + len_high_k
            len_high_k = length[j]
        fig, ax = plt.subplots()
        for i in range(1, nbd + 1):
            plt.plot(kpt, banddata[:, i], color=color_noispin)
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
        num_kpt, length, high_kpt_point, kpoints = get_kpt()
        nbd = get_NBANDS()
        banddata_up = np.loadtxt("band_up.dat", skiprows=1, dtype=np.float64)
        banddata_dw = np.loadtxt("band_dw.dat", skiprows=1, dtype=np.float64)
        kpt = banddata_up[:, 0]
        len_high_k = 0
        for j in range(0, len(length) - 1):
            length[j] = length[j] + len_high_k
            len_high_k = length[j]
        fig, ax = plt.subplots()
        for i in range(1, nbd + 1):
            plt.plot(kpt, banddata_up[:, i], color=color_noispin)
            plt.plot(kpt, banddata_dw[:, i], color=color_ispin)
        plt.xlim(0, max(length))
        ax.legend(("Spin_UP", "Spin_DW"), loc='best')
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
    else:
        print("No Ispin Information Obtained !")


# manipulate the fig
def manipulate_bandplot():
    print("    ****************************************************************")
    print("    *This is a code used to plot kinds of band structure,written by*")
    print("    *                     Xian-Yong Ding                           *")
    print("    ****************************************************************")
    print("\n")
    print("                       (^o^)GOOD LUCK!(^o^)                         ")
    print("\n")
    ispin = get_ispin()
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
                bandplot(energy, color)
            else:
                bandplot(energy, "black")
        elif inint == 2:
            color = str(input("which color do you want:"))
            print("(1) continue to setting energy range")
            print("(2) end setting")
            if con_end == 1:
                min_energy = float(input("minimum energy:"))
                max_energy = float(input("maximum energy:"))
                scale_energy = float(input("energy scale:"))
                energy = [min_energy, max_energy, scale_energy]
                bandplot(energy, color)
            else:
                bandplot([-10, 6, 2], color)
        else:
            bandplot([-10, 6, 2], "black")
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
                bandplot(energy, color_up, color_dw)
            else:
                bandplot(energy, "black", "red")
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
                bandplot(energy, color_up, color_dw)
            else:
                bandplot([-10, 6, 2], color_up, color_dw)
        else:
            bandplot([-10, 6, 2], "black", "red")
    else:
        print("    ******************** ISPIN is not 1 or 2!  Good bye ! ********************")



# *************************************** DOS *******************************************
# get element information from POSCAR
def _get_element():
    poscar = open('POSCAR', 'r')
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
def _dos_tot():
    len_zifu = 10
    len_zifu_cut = 12
    element, num_ele = _get_element()
    doscar = open('DOSCAR', 'r')
    dos = doscar.readlines()
    line_6 = dos[5].lstrip().split()
    eng_max = float(line_6[0])
    eng_min = float(line_6[1])
    nedos = int(line_6[2])
    fermi_eng = float(line_6[3])
    round(fermi_eng, 6)
    init = 1.00000
    ispin = get_ispin()
    if ispin == 1:
        with open("Tdos.dat", "w", encoding='utf-8') as tdos:
            tdos.write("Total density of state" + "\n")
            tdos.write("Energy".ljust(13," ") + "TDOS".ljust(13," ") + "IDOS".ljust(13," ") + "\n")
            for flag_line in range(0, nedos):
                tdos.write(str(round(float(dos[6+flag_line].split()[0])-fermi_eng, len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "   " +
                           str(round(float(dos[6+flag_line].split()[1]), len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "   "+
                           str(round(float(dos[6+flag_line].split()[2]), len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "\n")
        tdos.close()
    elif ispin == 2:
        with open("Tdos_SpinUp.dat", "w", encoding='utf-8') as tdos_up:
            tdos_up.write("Tdos_spinUP" + "\n")
            tdos_up.write("Energy".ljust(13," ") + "TDOS".ljust(13," ") + "IDOS".ljust(13," ") + "\n")
            for flag_line in range(0, nedos):
                tdos_up.write(str(round(float(dos[6+flag_line].split()[0])-fermi_eng, len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "   " +
                           str(round(float(dos[6+flag_line].split()[1]), len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "   "+
                           str(round(float(dos[6+flag_line].split()[3]), len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "\n")
        tdos_up.close()
        with open("Tdos_SpinDw.dat", "w", encoding='utf-8') as tdos_dw:
            tdos_dw.write("Tdos_spinDW" + "\n")
            tdos_dw.write("Energy".ljust(13," ") + "TDOS".ljust(13," ") + "IDOS".ljust(13," ") + "\n")
            for flag_line in range(0, nedos):
                tdos_dw.write(str(round(float(dos[6+flag_line].split()[0])-fermi_eng, len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "   " +
                           str(round(-float(dos[6+flag_line].split()[2]), len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "   "+
                           str(round(-float(dos[6+flag_line].split()[4]), len_zifu))[:len_zifu_cut].ljust(len_zifu_cut, " ") + "\n")
        tdos_dw.close()
    else:
        print("No ispin information obtained !")

# get pdos of atoms from DOSCAR
def _pdos_atom():
    orbitals = 10
    element, num_ele = _get_element()
    doscar = open('DOSCAR', 'r')
    dos = doscar.readlines()
    line_6 = dos[5].lstrip().split()
    eng_max = float(line_6[0])
    eng_min = float(line_6[1])
    nedos = int(line_6[2])
    fermi_eng = float(line_6[3])
    round(fermi_eng, 6)

    with open("tdos.dat", "w", encoding='utf-8') as tot:
        tot.write("Total density of state" + "\n")
        tot.write("energy  " + "TDOS   " + "     IDOS" + "\n")
        for flag_line in range(0, nedos):
            tot.write(str(dos[6 + flag_line].lstrip()))
    tot.close()
    file_dos = np.loadtxt("tdos.dat", skiprows=2, dtype=str)

    ispin = get_ispin()
    data_dos = []
    for lines in range(nedos + 7, len(dos), nedos+1):
        for lines_data in range(0, nedos):
            data_dos.append(dos[lines_data + lines])

    dos_atom_type = [[]] * len(num_ele)

    mid_dos_data = []
    atom_type = []
    for flag_mid_dos_data in range(0, sum(num_ele)):
        if flag_mid_dos_data == sum(num_ele):
            break
        else:
            mid_dos_data.append(data_dos[(flag_mid_dos_data * nedos):((flag_mid_dos_data + 1) * nedos)])
    for atom_type_num in range(0, len(num_ele)):
        atom_type.append(mid_dos_data[:num_ele[atom_type_num]])


    if ispin == 1:
        s = [[] for i in range(0, len(num_ele))]; s_m = [0] * nedos
        py = [[] for i in range(0, len(num_ele))]; py_m = [0] * nedos
        pz = [[] for i in range(0, len(num_ele))]; pz_m = [0] * nedos
        px = [[] for i in range(0, len(num_ele))]; px_m = [0] * nedos
        dxy = [[] for i in range(0, len(num_ele))]; dxy_m = [0] * nedos
        dyz = [[] for i in range(0, len(num_ele))]; dyz_m = [0] * nedos
        dz2 = [[] for i in range(0, len(num_ele))]; dz2_m = [0] * nedos
        dxz = [[] for i in range(0, len(num_ele))]; dxz_m = [0] * nedos
        dx2 = [[] for i in range(0, len(num_ele))]; dx2_m = [0] * nedos

        for flag_atom_type in range(0, len(num_ele)):
            for flag_atomnum in range(0, num_ele[flag_atom_type]):
                for flag_ne in range(0, nedos):
                    s[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[1]))
                    py[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[2]))
                    pz[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[3]))
                    px[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[4]))
                    dxy[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[5]))
                    dyz[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[6]))
                    dz2[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[7]))
                    dxz[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[8]))
                    dx2[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[9]))

        middle_s = [0] * nedos
        middle_py = [0] * nedos
        middle_pz = [0] * nedos
        middle_px = [0] * nedos
        middle_dxy = [0] * nedos
        middle_dyz = [0] * nedos
        middle_dz2 = [0] * nedos
        middle_dxz = [0] * nedos
        middle_dx2 = [0] * nedos

        for flag_atom_type in range(0, len(num_ele)):
            for flag_atomnum in range(0, num_ele[flag_atom_type]):
                for flag_ne in range(0, nedos):
                    s_m[flag_ne] = (s[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_s[flag_ne])
                    middle_s[flag_ne] = s_m[flag_ne]
                    py_m[flag_ne] = (py[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_py[flag_ne])
                    middle_py[flag_ne]=(py_m[flag_ne])
                    pz_m[flag_ne] = (pz[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_pz[flag_ne])
                    middle_pz[flag_ne]=(pz_m[flag_ne])
                    px_m[flag_ne] = (px[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_px[flag_ne])
                    middle_px[flag_ne]=(px_m[flag_ne])
                    dxy_m[flag_ne] = (dxy[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dxy[flag_ne])
                    middle_dxy[flag_ne]=(dxy_m[flag_ne])
                    dyz_m[flag_ne] = (dyz[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dyz[flag_ne])
                    middle_dyz[flag_ne]=(dyz_m[flag_ne])
                    dz2_m[flag_ne] = (dz2[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dz2[flag_ne])
                    middle_dz2[flag_ne]=(dz2_m[flag_ne])
                    dxz_m[flag_ne] = (dxz[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dxz[flag_ne])
                    middle_dxz[flag_ne]=(dxz_m[flag_ne])
                    dx2_m[flag_ne] = (dx2[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dx2[flag_ne])
                    middle_dx2[flag_ne]=(dx2_m[flag_ne])
        for ele in range(0, len(element)):
            with open("pdos_" + str(element[ele]) + ".dat", "w", encoding='utf-8') as pdos:
                pdos.write("pdos_" + str(element[ele]) +  "\n")
                pdos.write("Energy  s  py  pz  px  dxy  dyz  dz2  dxz  dx2" + "\n")
                for flag_line in range(0, nedos):
                    file_dos[flag_line, 0] = float(file_dos[flag_line, 0]) - fermi_eng
                    round(float(file_dos[flag_line, 0]), 6)
                    round(s_m[flag_line], 8)
                    round(py_m[flag_line], 8)
                    round(pz_m[flag_line], 8)
                    round(px_m[flag_line], 8)
                    round(dxy_m[flag_line], 8)
                    round(dyz_m[flag_line], 8)
                    round(dz2_m[flag_line], 8)
                    round(dxz_m[flag_line], 8)
                    round(dx2_m[flag_line], 8)
                    pdos.write(str(file_dos[flag_line, 0]) + "   " + str(s_m[flag_line]) + "   " + str(
                        py_m[flag_line]) +  "   " + str(
                        pz_m[flag_line]) +  "   " + str(
                        px_m[flag_line]) +  "   " + str(
                        dxy_m[flag_line]) +  "   " + str(
                        dyz_m[flag_line]) +  "   " + str(
                        dz2_m[flag_line]) +  "   " + str(
                        dxz_m[flag_line]) +  "   " + str(
                        dx2_m[flag_line]) + "\n")
            pdos.close()
    elif ispin == 2:
        s_atom_up = [[] for i in range(0, len(num_ele))]; s_up = [0] * nedos
        s_atom_dw = [[] for i in range(0, len(num_ele))]; s_dw = [0] * nedos
        py_atom_up = [[] for i in range(0, len(num_ele))]; py_up = [0] * nedos
        py_atom_dw = [[] for i in range(0, len(num_ele))]; py_dw = [0] * nedos
        pz_atom_up = [[] for i in range(0, len(num_ele))]; pz_up = [0] * nedos
        pz_atom_dw = [[] for i in range(0, len(num_ele))]; pz_dw = [0] * nedos
        px_atom_up = [[] for i in range(0, len(num_ele))]; px_up = [0] * nedos
        px_atom_dw = [[] for i in range(0, len(num_ele))]; px_dw = [0] * nedos
        dxy_atom_up = [[] for i in range(0, len(num_ele))]; dxy_up = [0] * nedos
        dxy_atom_dw = [[] for i in range(0, len(num_ele))]; dxy_dw = [0] * nedos
        dyz_atom_up = [[] for i in range(0, len(num_ele))]; dyz_up = [0] * nedos
        dyz_atom_dw = [[] for i in range(0, len(num_ele))]; dyz_dw = [0] * nedos
        dz2_atom_up = [[] for i in range(0, len(num_ele))]; dz2_up = [0] * nedos
        dz2_atom_dw = [[] for i in range(0, len(num_ele))]; dz2_dw = [0] * nedos
        dxz_atom_up = [[] for i in range(0, len(num_ele))]; dxz_up = [0] * nedos
        dxz_atom_dw = [[] for i in range(0, len(num_ele))]; dxz_dw = [0] * nedos
        dx2_atom_up = [[] for i in range(0, len(num_ele))]; dx2_up = [0] * nedos
        dx2_atom_dw = [[] for i in range(0, len(num_ele))]; dx2_dw = [0] * nedos

        for flag_atom_type in range(0, len(num_ele)):
            for flag_atomnum in range(0, num_ele[flag_atom_type]):
                for flag_ne in range(0, nedos):
                    s_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[1]))
                    s_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[2]))

                    py_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[3]))
                    py_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[4]))
                    pz_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[5]))
                    pz_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[6]))
                    px_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[7]))
                    px_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[8]))

                    dxy_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[9]))
                    dxy_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[10]))
                    dyz_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[11]))
                    dyz_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[12]))
                    dz2_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[13]))
                    dz2_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[14]))
                    dxz_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[15]))
                    dxz_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[16]))
                    dx2_atom_up[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[17]))
                    dx2_atom_dw[flag_atom_type].append(float(atom_type[flag_atom_type][flag_atomnum][flag_ne].split()[18]))

        middle_s_up = [0] * nedos; middle_s_dw = [0] * nedos
        middle_py_up = [0] * nedos; middle_py_dw = [0] * nedos
        middle_pz_up = [0] * nedos; middle_pz_dw = [0] * nedos
        middle_px_up = [0] * nedos; middle_px_dw = [0] * nedos
        middle_dxy_up = [0] * nedos; middle_dxy_dw = [0] * nedos
        middle_dyz_up = [0] * nedos; middle_dyz_dw = [0] * nedos
        middle_dz2_up = [0] * nedos; middle_dz2_dw = [0] * nedos
        middle_dxz_up = [0] * nedos; middle_dxz_dw = [0] * nedos
        middle_dx2_up = [0] * nedos; middle_dx2_dw = [0] * nedos

        for flag_atom_type in range(0, len(num_ele)):
            for flag_atomnum in range(0, num_ele[flag_atom_type]):
                for flag_ne in range(0, nedos):
                    s_up[flag_ne] = (s_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_s_up[flag_ne])
                    middle_s_up[flag_ne] = s_up[flag_ne]
                    s_dw[flag_ne] = (s_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_s_dw[flag_ne])
                    middle_s_dw[flag_ne] = s_dw[flag_ne]

                    py_up[flag_ne] = (py_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_py_up[flag_ne])
                    middle_py_up[flag_ne]=(py_up[flag_ne])
                    py_dw[flag_ne] = (py_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_py_dw[flag_ne])
                    middle_py_dw[flag_ne] = (py_dw[flag_ne])

                    pz_up[flag_ne] = (pz_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_pz_up[flag_ne])
                    middle_pz_up[flag_ne]=(pz_up[flag_ne])
                    pz_dw[flag_ne] = (pz_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_pz_dw[flag_ne])
                    middle_pz_dw[flag_ne] = (pz_dw[flag_ne])

                    px_up[flag_ne] = (px_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_px_up[flag_ne])
                    middle_px_up[flag_ne]=(px_up[flag_ne])
                    px_dw[flag_ne] = (px_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_px_dw[flag_ne])
                    middle_px_dw[flag_ne] = (px_dw[flag_ne])

                    dxy_up[flag_ne] = (dxy_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dxy_up[flag_ne])
                    middle_dxy_up[flag_ne]=(dxy_up[flag_ne])
                    dxy_dw[flag_ne] = (dxy_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dxy_dw[flag_ne])
                    middle_dxy_dw[flag_ne] = (dxy_dw[flag_ne])

                    dyz_up[flag_ne] = (dyz_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dyz_up[flag_ne])
                    middle_dyz_up[flag_ne]=(dyz_up[flag_ne])
                    dyz_dw[flag_ne] = (dyz_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dyz_dw[flag_ne])
                    middle_dyz_dw[flag_ne] = (dyz_dw[flag_ne])

                    dz2_up[flag_ne] = (dz2_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dz2_up[flag_ne])
                    middle_dz2_up[flag_ne]=(dz2_up[flag_ne])
                    dz2_dw[flag_ne] = (dz2_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dz2_dw[flag_ne])
                    middle_dz2_dw[flag_ne] = (dz2_dw[flag_ne])

                    dxz_up[flag_ne] = (dxz_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dxz_up[flag_ne])
                    middle_dxz_up[flag_ne]=(dxz_up[flag_ne])
                    dxz_dw[flag_ne] = (dxz_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dxz_dw[flag_ne])
                    middle_dxz_dw[flag_ne] = (dxz_dw[flag_ne])

                    dx2_up[flag_ne] = (dx2_atom_up[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dx2_up[flag_ne])
                    middle_dx2_up[flag_ne]=(dx2_up[flag_ne])
                    dx2_dw[flag_ne] = (dx2_atom_dw[flag_atom_type][flag_ne + flag_atomnum * nedos] + middle_dx2_dw[flag_ne])
                    middle_dx2_dw[flag_ne] = (dx2_dw[flag_ne])

        for ele in range(0, len(element)):
            with open("pdos_" + str(element[ele]) + ".dat", "w", encoding='utf-8') as pdos:
                pdos.write("pdos_" + str(element[ele]) +  "\n")
                pdos.write("Energy  s_up  s_dw  py_up  py_dw  pz_up  pz_dw  px_up  px_dw  dxy_up  dxy_dw  dyz_up  dyz_dw  dz2_up  dz2_dw  dxz_up  dxz_dw  dx2_up  dx2_dw" + "\n")
                for flag_line in range(0, nedos):
                    file_dos[flag_line, 0] = float(file_dos[flag_line, 0]) - fermi_eng
                    round(float(file_dos[flag_line, 0]), 6)

                    round(s_up[flag_line], 8)
                    round(s_dw[flag_line], 8)

                    round(py_up[flag_line], 8)
                    round(py_dw[flag_line], 8)
                    round(pz_up[flag_line], 8)
                    round(pz_dw[flag_line], 8)
                    round(px_up[flag_line], 8)
                    round(px_dw[flag_line], 8)

                    round(dxy_up[flag_line], 8)
                    round(dxy_dw[flag_line], 8)
                    round(dyz_up[flag_line], 8)
                    round(dyz_dw[flag_line], 8)
                    round(dz2_up[flag_line], 8)
                    round(dz2_dw[flag_line], 8)
                    round(dxz_up[flag_line], 8)
                    round(dxz_dw[flag_line], 8)
                    round(dx2_up[flag_line], 8)
                    round(dx2_dw[flag_line], 8)
                    pdos.write(str(file_dos[flag_line, 0]) + "   " + str(
                        s_up[flag_line]) + "   " + str(
                        -s_dw[flag_line]) + "   " + str(
                        py_up[flag_line]) + "   " + str(
                        -py_dw[flag_line]) + "   " + str(
                        pz_up[flag_line]) + "   " + str(
                        -pz_dw[flag_line]) + "   " + str(
                        px_up[flag_line]) + "   " + str(
                        -px_dw[flag_line]) + "   " + str(
                        dxy_up[flag_line]) + "   " + str(
                        -dxy_dw[flag_line]) + "   " + str(
                        dyz_up[flag_line]) + "   " + str(
                        -dyz_dw[flag_line]) + "   " + str(
                        dz2_up[flag_line]) + "   " + str(
                        -dz2_dw[flag_line]) + "   " + str(
                        dxz_up[flag_line]) + "   " + str(
                        -dxz_dw[flag_line]) + "   " + str(
                        dx2_up[flag_line]) + "   " + str(
                        -dx2_dw[flag_line]) +
                                  "\n")
            pdos.close()
    else:
        print("No Ispin Information Obtained !, Checking whether this path exist the INCAR file")
# dos plot
def _dos_plot_atom_orbital(eng_range = [-10, 10, 2]):
    _pdos_atom()
    orbitals_noispin = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"]
    orbitals_ispin = ["s_up", "s_dw", "py_up", "py_dw", "pz_up", "pz_dw", "px_up", "px_dw", "dxy_up", "dxy_dw", "dyz_up", "dyz_dw", "dz2_up", "dz2_dw", "dxz_up", "dxz_dw", "dx2_up", "dx2_dw"]
    ispin = get_ispin()
    element, num_ele = _get_element()
    filename = [[] for i in range(0, len(element))]
    data = [[] for i in range(0, len(element))]
    for flag_atom_type in range(0, len(element)):
        filename[flag_atom_type] = "pdos_" + element[flag_atom_type] + ".dat"
        data[flag_atom_type] = np.loadtxt(filename[flag_atom_type], skiprows=2, dtype=float)

    get_d_maxvalue = [[] for i in range(0, len(element))]
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    if ispin == 1:
        print(" *************************** ISPIN is 1 ****************************")
        for flag_atom_type in range(0, len(element)):
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 6])
            if get_d_maxvalue[flag_atom_type] == 0:     # if d orbitals have no values, plot the s p
                for flag_atom_orbital in range(0, (len(orbitals_noispin)-5)):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, 1+flag_atom_orbital], ls="-",
                             label=str(element[flag_atom_type] +"_"+ orbitals_noispin[flag_atom_orbital]))
            else:
                for flag_atom_orbital in range(0, len(orbitals_noispin)):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, 1+flag_atom_orbital], ls="-",
                             label=str(element[flag_atom_type] +"_"+ orbitals_noispin[flag_atom_orbital]))

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
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 10])
            if get_d_maxvalue[flag_atom_type] == 0:     # if d orbitals have no values, plot the s p
                for flag_atom_orbital in range(0, (len(orbitals_ispin)-10)):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, (1+flag_atom_orbital)], ls="-",
                             label=str(element[flag_atom_type] +"_"+ orbitals_ispin[flag_atom_orbital]))
            else:
                for flag_atom_orbital in range(0, len(orbitals_ispin)):
                    plt.plot(data[flag_atom_type][:, 0], data[flag_atom_type][:, 1 + flag_atom_orbital], ls="-",
                             label=str(element[flag_atom_type] + "_" + orbitals_ispin[flag_atom_orbital]))

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



# plot dos for atom_orbital
def _dos_plot_atom(eng_range = [-10, 10, 2]):
    _pdos_atom()
    orbitals_noispin = ["s", "p", "d"]
    orbitals_ispin = ["s_up", "s_dw", "p_up", "p_dw", "d_up", "d_dw"]
    ispin = get_ispin()
    element, num_ele = _get_element()
    filename = [[] for i in range(0, len(element))]
    data = [[] for i in range(0, len(element))]
    for flag_atom_type in range(0, len(element)):
        filename[flag_atom_type] = "pdos_" + element[flag_atom_type] + ".dat"
        data[flag_atom_type] = np.loadtxt(filename[flag_atom_type], skiprows=2, dtype=float)

    get_d_maxvalue = [[] for i in range(0, len(element))]
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    if ispin == 1:
        print("    *************************** ISPIN is 1 ****************************")
        s_orbital = [[] for i in range(0, len(element))]
        p_orbital = [[] for i in range(0, len(element))]
        d_orbital = [[] for i in range(0, len(element))]
        for flag_atomtype in range(0, len(element)):
            s_orbital[flag_atomtype] = data[flag_atomtype][:, 1]
            p_orbital[flag_atomtype] = data[flag_atomtype][:, 2] + data[flag_atomtype][:, 3] + data[flag_atomtype][:, 4]
            d_orbital[flag_atomtype] = data[flag_atomtype][:, 5] + data[flag_atomtype][:, 6] + data[flag_atomtype][:, 7] + data[flag_atomtype][:, 8] + data[flag_atomtype][:, 9]

        for flag_atom_type in range(0, len(element)):
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 6])
            if get_d_maxvalue[flag_atom_type] == 0:     # if d orbitals have no values, plot the s p
                plt.plot(data[flag_atom_type][:, 0], s_orbital[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_noispin[0]))
                plt.plot(data[flag_atom_type][:, 0], p_orbital[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_noispin[1]))
            else:
                plt.plot(data[flag_atom_type][:, 0], s_orbital[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_noispin[0]))
                plt.plot(data[flag_atom_type][:, 0], p_orbital[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_noispin[1]))
                plt.plot(data[flag_atom_type][:, 0], d_orbital[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_noispin[2]))

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
        s_orbital_up = [[] for i in range(0, len(element))]
        s_orbital_dw = [[] for i in range(0, len(element))]
        p_orbital_up = [[] for i in range(0, len(element))]
        p_orbital_dw = [[] for i in range(0, len(element))]
        d_orbital_up = [[] for i in range(0, len(element))]
        d_orbital_dw = [[] for i in range(0, len(element))]
        for flag_atomtype in range(0, len(element)):
            s_orbital_up[flag_atomtype] = data[flag_atomtype][:, 1]
            s_orbital_dw[flag_atomtype] = data[flag_atomtype][:, 2]

            p_orbital_up[flag_atomtype] = data[flag_atomtype][:, 3] + data[flag_atomtype][:, 5] + data[flag_atomtype][:, 7]
            p_orbital_dw[flag_atomtype] = data[flag_atomtype][:, 4] + data[flag_atomtype][:, 6] + data[flag_atomtype][:, 8]

            d_orbital_up[flag_atomtype] = data[flag_atomtype][:, 9] + data[flag_atomtype][:, 11] + data[flag_atomtype][:, 13] + \
                                       data[flag_atomtype][:, 15] + data[flag_atomtype][:, 17]
            d_orbital_dw[flag_atomtype] = data[flag_atomtype][:, 10] + data[flag_atomtype][:, 12] + data[flag_atomtype][:, 14] + \
                                       data[flag_atomtype][:, 16] + data[flag_atomtype][:, 18]
        for flag_atom_type in range(0, len(element)):
            get_d_maxvalue[flag_atom_type] = max(data[flag_atom_type][:, 6])
            if get_d_maxvalue[flag_atom_type] == 0:  # if d orbitals have no values, plot the s p
                plt.plot(data[flag_atom_type][:, 0], s_orbital_up[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[0]))
                plt.plot(data[flag_atom_type][:, 0], s_orbital_dw[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[1]))

                plt.plot(data[flag_atom_type][:, 0], p_orbital_up[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[2]))
                plt.plot(data[flag_atom_type][:, 0], p_orbital_dw[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[3]))

            else:
                plt.plot(data[flag_atom_type][:, 0], s_orbital_up[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[0]))
                plt.plot(data[flag_atom_type][:, 0], s_orbital_dw[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[1]))

                plt.plot(data[flag_atom_type][:, 0], p_orbital_up[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[2]))
                plt.plot(data[flag_atom_type][:, 0], p_orbital_dw[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[3]))

                plt.plot(data[flag_atom_type][:, 0], d_orbital_up[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[4]))
                plt.plot(data[flag_atom_type][:, 0], d_orbital_dw[flag_atom_type], ls="-",
                         label=str(element[flag_atom_type] + "_" + orbitals_ispin[5]))

        ax.axvline(x=0, ymin=0, ymax=1, linestyle='--', linewidth=1.5, color='0.5')
        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1.5, color='0.5')
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


# plot tot dos
def _dos_plot_tot(eng_range = [-10, 10, 2], integral = 0):
    _dos_tot()
    ispin = get_ispin()
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

# manipulate dos ploting
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


# ***************************** main control ************************************
def main():
    print("    ********************************************************************")
    print("    *This is a code used to simplify the vasp calculation, written by **")
    print("    ******************           XY Ding                ****************")
    print("    ********************************************************************")
    print("    *************          (^o^)GOOD LUCK!(^o^)           **************")
    print("\n")
    print("    *************** post-processing for vasp calculation ***************")
    print("********  Data processing for vasp calculation: band and dos ***********")
    print("(1) band " + "\t" + "(2) dos ")
    process_num = int(input("Input a Number: "))
    if process_num == 1:
        manipulate_bandplot()
    elif process_num == 2:
        _manipulate_dos()
    else:
        print("You are input a wrong number")

if __name__ == '__main__':
    main()



