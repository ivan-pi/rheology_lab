import csv
import numpy as np

def podatki(filename,stVrst):

    # INITIALIZE LISTS
    col1 = []
    col2 = []
    col3 = []
    col4 = []

    # READ FROM CSV FILE
    with open(filename, 'rb') as csvfile:
        for i in range(0,stVrst):
            next(csvfile)
        podatki = csv.reader(csvfile)
        for row in podatki:
            col1.append(float(row[1]))
            col2.append(float(row[2]))
            col3.append(float(row[3]))
            col4.append(float(row[4]))

    # CONVERT LIST TO ARRAYS
    col1 = np.array(col1)
    col2 = np.array(col2)
    col3 = np.array(col3)
    col4 = np.array(col4)

    return col1,col2,col3,col4


# PROGRAM STARTS HERE

file1 = "S1TempX1.0.csv"
file2 = "S2TempG1.0.csv"

temp1,stress1,visc1,speed1  = podatki(file1,6)
shear2,temp2,stress2,visc2 = podatki(file2,2)

print len(temp2)
#-------------------------------
# GRAF
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import gridspec

red_patch = mpatches.Patch(color='red', label='gelan 0.7 %, 20 $^\circ$C')
blue_patch = mpatches.Patch(color='blue', label='gelan 1.0 %, 20 $^\circ$C')
green_patch = mpatches.Patch(color='green', label='xantan 1.0 %, 20 $^\circ$C')

T_dots = mlines.Line2D([], [], color='blue', linestyle='-',marker='o',label="T $\uparrow$")
T_diam = mlines.Line2D([], [], color='blue', linestyle='--',marker='d',label="T $\downarrow$")

plt.figure(figsize=(8,4))
plt.gcf().subplots_adjust(bottom=0.15)

# TEST LEZENJA IN OBNOVE
plt.plot(temp2[0:30],visc2[0:30],'bo-',
         temp2[30:60],visc2[30:60],'bd-')

plt.legend(handles=[T_dots,T_diam],loc='upper right',fontsize='14',ncol=1)
plt.title('gelan 1,0 %, $\dot{\\gamma}$ = 50 s$^{-1}$',fontsize='16')
plt.ylabel('$\eta$ (Pa$\cdot$s)',fontsize='16')
plt.xlabel('T ($^\circ$C)', fontsize='16')
plt.yscale("log")

plt.savefig("temperature.png")
plt.show()
plt.close()
