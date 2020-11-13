# -*- coding: utf-8 -*-
"""
Created on 10:40 12-11-2020 

@author: Xian-Yong Ding
mail to: dxy_vasp@163.com
python3: optics.py
"""
import numpy as np
from matplotlib import pyplot as plt
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import math 

# get NEDOS from OUTCAR file
def get_nedos(filename = "OUTCAR"):
    outcar = open(filename, 'r')
    outcar_lines = outcar.readlines()
    for flag_lines in range(0, len(outcar_lines)):
        if 'NEDOS' in outcar_lines[flag_lines]:
            nedos = outcar_lines[flag_lines].split('=')[1].split()[0]
    return int(nedos)

# obtain data of real and imaginary part
def get_data(filename_real = "real.dat", filename_imag = "imag.dat"):
    nedos = get_nedos()
    h = 4.1356676969 * (10 ** (-15)) / (2 * np.pi)
    data_real = open(filename_real, 'r')
    data_imag = open(filename_imag, 'r')
    data_real_lines = data_real.readlines()
    data_imag_lines = data_imag.readlines()
    energy = np.zeros(shape=(nedos))
    w1 = np.zeros(shape=(nedos))
    w2 = np.zeros(shape=(nedos))
    for flag_lines in range(0, nedos):
        energy[flag_lines] = data_imag_lines[flag_lines].split()[0]
        w1[flag_lines] = float(data_real_lines[flag_lines].split()[1])
        w2[flag_lines] = float(data_imag_lines[flag_lines].split()[1])
    n_w = np.zeros(shape=(nedos))
    k_w = np.zeros(shape=(nedos))
    a_w = np.zeros(shape=(nedos))
    l_w = np.zeros(shape=(nedos))
    r_w = np.zeros(shape=(nedos))
    for i in range(0, nedos):
        n_w[i] = ((((w1[i] ** 2 + w2[i] ** 2)**(1/2) + w1[i])/2) ** (1/2))
        omega = energy[i]/h
        k_w[i] = ((((w1[i] ** 2 + w2[i] ** 2)**(1/2) - w1[i])/2) ** (1/2))
        a_w[i] = omega * (2 ** (1/2)) * ((((w1[i] ** 2 + w2[i] ** 2)**(1/2) - w1[i])) ** (1/2))/(3*(10**8) * 100)
        l_w[i] = w2[i]/(w1[i] ** 2 + w2[i] ** 2)
        r_w[i] = ((n_w[i]-1) ** 2 + k_w[i] ** 2)/((n_w[i] + 1) ** 2 + k_w[i] ** 2)
    return energy, n_w, k_w, a_w, l_w, r_w

# plot
def plot_data(x, y, ylabel, color = "blue", energy_gr = [0, 24, 2]):
    min_value, max_value = get_max(x, y)
    plt.plot(x, y, color=color)
    plt.xlim(energy_gr[0], energy_gr[1])
    plt.ylim(min_value, max_value*1.1)
    plt.xlabel("Energy (eV)", fontsize=14, fontname='arial')
    plt.ylabel(ylabel, fontsize=14, fontname='arial')
    xtick = np.arange(energy_gr[0], energy_gr[1] + 1, energy_gr[2])
    a = int(len(xtick) / 2)
    plt.xticks(np.insert(xtick, a, 0), fontsize=12)
    plt.savefig(ylabel+".png", dpi=300)
    plt.show()
def set_plot(x, y, type_plot):
    print("    ************************** Entering" + str(type_plot) + "index plot ****************************")
    print("(1) Setting energy")
    print("(2) Setting color")
    print("(3) Use default setting")
    select = int(input("please input a number: "))
    if select == 1:
        min_energy = float(input("minimum energy:"))
        max_energy = float(input("maximum energy:"))
        scale_energy = float(input("energy scale:"))
        energy_gr = [min_energy, max_energy, scale_energy]
        print("(1) Continue to setting color")
        print("(2) End setting")
        select_go = int(input("please input a number: "))
        if select_go == 1:
            color = str(input("color is: "))
            plot_data(x,y, type_plot, color, energy_gr)
        else:
            plot_data(x, y, type_plot, "blue", energy_gr)
    elif select == 2:
        color = str(input("color is: "))
        print("(1) Continue to setting energy")
        print("(2) End setting")
        select_go = int(input("please input a number: "))
        if select_go == 1:
            min_energy = float(input("minimum energy:"))
            max_energy = float(input("maximum energy:"))
            scale_energy = float(input("energy scale:"))
            energy_gr = [min_energy, max_energy, scale_energy]
        else:
            plot_data(x, y, type_plot, color)
    else:
        plot_data(x, y, type_plot)
# plot all data in one figure
def get_max(x_range, col_type):
    data = []
    for i in range(0, len(x_range)):
        if 0 <= x_range[i] <= 24:
            data.append(col_type[i])
    max_value = max(data)
    min_value = min(data)- min(data) * 0.1
    return min_value, max_value

def set_plot_all(energy, colume, color_type = ["red", "blue", "green", "black", "orange"], energy_gr=[0, 24, 4]):
    type_plot=["Refractive", "Extinction", "Absorption (cm$^{-1}$)", "Energy loss", "Reflectivity"]
    plt.figure(figsize=(20, 16), dpi=80)
    for i in range(1, 6):
        min_value, max_value = get_max(energy, colume[i-1])
        ax = plt.subplot(3, 2, i)
        ax.plot(energy, colume[i-1], color=color_type[i-1])
        ax.set_xlim(energy_gr[0], energy_gr[1])
        ax.set_ylim(min_value, max_value*1.1)
        ax.set_xlabel("Energy (eV)", fontsize=14, fontname='arial')
        ax.set_ylabel(type_plot[i-1], fontsize=14, fontname='arial')
        xtick = np.arange(energy_gr[0], energy_gr[1] + 1, energy_gr[2])
        a = int(len(xtick) / 2)
        ax.set_xticks(np.insert(xtick, a, 0))
    plt.savefig("All.png", dpi=300)
    plt.show()
# manipulate plot
def manipulate_plot():
    energy, refractive, extinction, absorption, energy_loss, reflectivity = get_data()
    print("    ************************** which type of plot do you want ? ****************************")
    print("(1) Refractive index")
    print("(2) Extinction coefficient")
    print("(3) Absorption coefficient")
    print("(4) Energy loss function")
    print("(5) reflectivity")
    print("(6) All of them")
    select_type = int(input("please input a number: "))
    if select_type == 1:
        set_plot(energy, refractive, "Refractive")
    elif select_type == 2:
        set_plot(energy, extinction, "Extinction")
    elif select_type == 3:
        set_plot(energy, absorption, "Absorption (cm$^{-1}$)")
    elif select_type == 4:
        set_plot(energy, energy_loss, "Energy_loss")
    elif select_type == 5:
        set_plot(energy, reflectivity, "Reflectivity")
    elif select_type == 6:
        colume = [refractive, extinction, absorption, energy_loss, reflectivity]
        set_plot_all(energy, colume)
    else:
        print("you are input a wrong number, please re-input !")

if __name__ == '__main__':
    manipulate_plot()