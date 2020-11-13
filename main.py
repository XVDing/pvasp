# -*- coding: utf-8 -*-
"""
Created on 10:21 14-10-2020 

@author: Xian-Yong Ding
mail to: dxy_vasp@163.com
python3: main.py
"""
import numpy as np
from matplotlib import pyplot as plt
import math
# import other python scripts
import os, sys
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)

import scripts.initial
import scripts.band as bd
import scripts.hse06 as hse
import scripts.dos as dos
import scripts.optics as op


file_path = curPath + "/INCAR_templates"

# main control
print(" ********************************************************************")
print(" *This is a code used to simplify the vasp calculation, written by **")
print(" *************           Xian-Yong Ding                **************")
print(" Thanks for your citing: New Journal of Physics, Volume 22, July 2020")
print(" ********************************************************************")
print(" *************          (^o^)GOOD LUCK!(^o^)           **************")
print("\n")
print(" *********  Pre and post-preparation for vasp calculation  **********")
print(" **********  eg:(1) pre-processing: incar template    ***************")
print(" **********  eg:(2) post-processing: band plot    *******************")
print(" ********************************************************************")
print(" (1) Incar template for vasp calculation")
print(" (2) Data processing for vasp calculation:eg. band and dos")
print(" (3) User defined plot scripts, data from first line ")
seltct_vasp = int(input("Input a number: "))
print("\n")
if seltct_vasp == 1:
    print("********  Incar template for vasp calculation !  *******")
    print("(1) optimization" + "\t" + "(2) self-consistent")
    print("(3) band calculation" + "\t" + "(4) dos calculation")
    print("(5) HSE06 calculation" + "\t" + "(6) phonon")
    print("(7) Elastic calculation" + "\t" + "(8) md")
    incar = int(input("Input a Number: "))
    if incar == 1:
        os.system("cp " + file_path +"/INCAR_opt ./INCAR")
    elif incar == 2:
        os.system("cp " + file_path +"/INCAR_scf ./INCAR")
    elif incar == 3:
        os.system("cp " + file_path +"/INCAR_band ./INCAR")
    elif incar == 4:
        os.system("cp " + file_path +"/INCAR_dos ./INCAR")
    elif incar == 5:
        os.system("cp " + file_path +"/INCAR_hse ./INCAR")
    elif incar == 6:
        os.system("cp " + file_path +"/INCAR_phonon ./INCAR")
    elif incar == 7:
        os.system("cp " + file_path +"/INCAR_elastic ./INCAR")
    elif incar == 8:
        os.system("cp " + file_path +"/INCAR_md ./INCAR")
    else:
        print("you have print a wrong number !")
elif seltct_vasp == 2:
    print("********  Data processing for vasp calculation: band and dos ********")
    print("(1) band plot       " + "\t" + "(2) HSE06 band plot ")
    print("(3) Projected band  " + "\t" + "(4) Density of state ")
    print("(5) Optics properties")
    process_num = int(input("Input a Number: "))
    if process_num == 1:
        print(" ************************** Enter band ploting ****************************")
        bd.manipulate_bandplot()
    elif process_num == 2:
        print(" ************************** Enter HSE06 ploting ****************************")
        hse.manipulate_bandplot()
    elif process_num == 3:
        print(" ****************************** waitting ***********************************")
    elif process_num == 4:
        print(" ************************** Enter DOS ploting ******************************")
        dos._manipulate_dos()
    elif process_num == 5:
        print(" ******************** Enter optics properties ploting **********************")
        os.system("cp " + curPath + "/scripts/optics.sh ./")
        os.system("bash optics.sh")
        op.manipulate_plot()
    else:
        print("please input a right number !")
elif seltct_vasp == 3:
    print("********  Data processing for User defined scripts: data from first line, two colume ********")
else:
    print("You have input a wrong number !")




