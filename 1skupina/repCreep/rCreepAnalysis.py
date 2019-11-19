import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import rCreep as rcr

# imena datotek s podatki
file1 = 'S1gelan0.7.csv'
file2 = 'S1gelan1.0.csv'
file3 = 'S1ksantan1.0.csv'

file5 = 'S2gelan1.0.csv'
file6 = 'S2ksantan1.0.csv'

# uvozimo podatke
time1, shear1, strain1, comp1 = rcr.podatki(file1, 4)
time2, shear2, strain2, comp2 = rcr.podatki(file2, 5)
time3, shear3, strain3, comp3 = rcr.podatki(file3, 4)

time5, shear5, strain5, comp5 = rcr.podatki(file5, 4)
time6, shear6, strain6, comp6 = rcr.podatki(file6, 6)

print len(strain1), len(strain2), len(strain3), len(strain5), len(strain6)

tauC = 0.03
Jnr = 0.0
R = 0.0
gnr = 0
gp = 0
for i in range(0, 3):
    gp = strain1[i*30+9] - gnr
    gnr = strain1[i*30+29] - gnr

    Jnr = Jnr + gnr
    R = R + (gp-gnr)/gp

Jnr = Jnr/(tauC*4.0)
R = R/4 * 100

print "gelan 0.7 1. skup,","Jnr =",Jnr,",R =", R

tauC = 0.3
Jnr = 0.0
R = 0.0
gnr = 0
gp = 0
for i in range(0, 3):
    gp = strain2[i*30+9] - gnr
    gnr = strain2[i*30+29] - gnr

    Jnr = Jnr + gnr
    R = R + (gp-gnr)/gp

Jnr = Jnr/(tauC*4.0)
R = R/4 * 100

print "gelan 1,0 1. skup","Jnr =",Jnr,",R =", R

tauC = 0.1
Jnr = 0.0
R = 0.0
gnr = 0
gp = 0
for i in range(0, 1):
    gp = strain3[i*30+9] - gnr
    gnr = strain3[i*30+29] - gnr

    Jnr = Jnr + gnr
    R = R + (gp-gnr)/gp

Jnr = Jnr/(tauC*2.0)
R = R/2 * 100

print "ksantan 1,0 1. skup","Jnr =",Jnr,",R =", R

#---------------------------------------------

tauC = 0.3
Jnr = 0.0
R = 0.0
gnr = 0
gp = 0
for i in range(0, 3):
    gp = strain5[i*30+9] - gnr
    gnr = strain5[i*30+29] - gnr

    Jnr = Jnr + gnr
    R = R + (gp-gnr)/gp

Jnr = Jnr/(tauC*4.0)
R = R/4 * 100

print "gelan 1,0 2. skup","Jnr =",Jnr,",R =", R

tauC = 0.1
Jnr = 0.0
R = 0.0
gnr = 0
gp = 0
for i in range(0, 3):
    gp = strain6[i*30+9] - gnr
    gnr = strain6[i*30+29] - gnr

    Jnr = Jnr + gnr
    R = R + (gp-gnr)/gp

Jnr = Jnr/(tauC*4.0)
R = R/4 * 100

print "ksantan 1,0 2. skup","Jnr =",Jnr,",R =", R

#-------------------------------
# GRAF
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import gridspec

red_patch = mpatches.Patch(color='red', label='gelan 0.7 %, 20 $^\circ$C')
blue_patch = mpatches.Patch(color='blue', label='gelan 1.0 %, 20 $^\circ$C')
green_patch = mpatches.Patch(color='green', label='xantan 1.0 %, 20 $^\circ$C')

J_dots = mlines.Line2D([], [], color='grey', linestyle='-',marker='o',label="$\\gamma$, 1.skupina")
J_diam = mlines.Line2D([], [], color='grey', linestyle='--',marker='d',label="$\\gamma$, 2.skupina")

plt.figure(figsize=(15,6))
plt.gcf().subplots_adjust(bottom=0.15)
plt.subplots_adjust(wspace=0.3)

# TEST LEZENJA IN OBNOVE
plt.subplot(121)
plt.plot(time1,strain1,'ro-',
         time2,strain2,'bo-',
         time5,strain6,'bd--')

plt.legend(handles=[red_patch,blue_patch,J_dots,J_diam],loc='upper left',fontsize='11',ncol=1)
plt.title('gelan',fontsize='18')
plt.ylabel('$\\gamma$ (/)',fontsize='16')
plt.xlabel('t (s)', fontsize='16')

plt.subplot(122)
plt.plot(time3,strain3,'go-',
         time6,strain6,'gd--')

plt.legend(handles=[green_patch,J_dots,J_diam],loc='lower right',fontsize='11',ncol=1)
plt.title('ksantan',fontsize='18')
plt.ylabel('$\\gamma$ (/)',fontsize='16')
plt.xlabel('t (s)', fontsize='16')

plt.savefig("rCreep.png")
plt.show()
plt.close()
