import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import amplTest as at

# stevilo praznih vrstic
prazneVrstice = 5

# imena datotek s podatki
file1 = 'aS2gelanP1.0T30.csv'

file4 = 'aS2xantanP1T20.csv'
file5 = 'aS2xantanP1T30.csv'

# uvozimo podatke
strain1, shear1, G1_1, G2_1 = at.podatki(file1, prazneVrstice)

strain4, shear4, G1_4, G2_4 = at.podatki(file4, prazneVrstice)
strain5, shear5, G1_5, G2_5 = at.podatki(file5, prazneVrstice)

# izracunamo kompleksni modul
Gc1 = np.sqrt(G1_1**2+G2_1**2)

Gc4 = np.sqrt(G1_4**2+G2_4**2)
Gc5 = np.sqrt(G1_5**2+G2_5**2)

# nastavimo zacetne vrednosti in omejitve parametrov
par1 = Parameters()
par1.add('G0',value=1,min=0.0)
par1.add('a' ,value=1,min=0.0)
par1.add('b' ,value=1,min=0.0)
par1.add('n' ,value=1,min=0.0)

# minimizacija napake - non-linear least squares method
out1 = minimize(at.amplModel,par1,args=(strain1,Gc1),method='leastsq')

# pridobimo parametre iz fitanja
k1 = at.koeficienti(par1)

# nastavimo zacetne vrednosti in omejitve parametrov
par4 = Parameters()
par4.add('G0',value=1,min=0.0)
par4.add('a' ,value=1,min=0.0)
par4.add('b' ,value=1,min=0.0)
par4.add('n' ,value=1,min=0.0)

# minimizacija napake - non-linear least squares method
out4 = minimize(at.amplModel,par4,args=(strain4,Gc4),method='leastsq')

# pridobimo parametre iz fitanja
k4 = at.koeficienti(par4)

# nastavimo zacetne vrednosti in omejitve parametrov
par5 = Parameters()
par5.add('G0',value=1,min=0.0)
par5.add('a' ,value=1,min=0.0)
par5.add('b' ,value=1,min=0.0)
par5.add('n' ,value=1,min=0.0,max=5.0)

# minimizacija napake - non-linear least squares method
out5 = minimize(at.amplModel,par5,args=(strain5,Gc5),method='leastsq')

# pridobimo parametre iz fitanja
k5 = at.koeficienti(par5)

krit1 = (0.03/(0.97*k1[1]-k1[2]))**(1.0/k1[3])
krit4 = (0.03/(0.97*k4[1]-k4[2]))**(1.0/k4[3])
krit5 = (0.03/(0.97*k5[1]-k5[2]))**(1.0/k5[3])

# natisnemo parametre na zaslon
print "k1",krit1, "gelan 1.0, 30C", ["%0.12f" % i for i in k1]
print "k4",krit4, "xantan 20C", ["%0.12f" % i for i in k4]
print "k5",krit5, "xantan 30C", ["%0.12f" % i for i in k5]

# GRAF
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import rc,rcParams

red_patch   = mpatches.Patch(color='red',  label='gelan 1 %, 30 $^\circ$C')
red_patch2  = mpatches.Patch(color='red',  label='ksantan 1 %, 20 $^\circ$C')
blue_patch2 = mpatches.Patch(color='blue', label='ksantan 1 %, 30 $^\circ$C')
g1_line = mlines.Line2D([], [], color='grey', marker='o',label="G*")
g2_dots = mlines.Line2D([], [], color='grey', marker='d',label="$\\tau$")
krit_line = mlines.Line2D([], [], color='grey', linestyle='--',label="$\\gamma_{\mathrm{krit}}$")

plt.figure(figsize=(15,5))
plt.gcf().subplots_adjust(bottom=0.15)
plt.subplots_adjust(wspace=0.4)

# MEHANSKI SPEKTER
ax1 = plt.subplot(121)

x = np.logspace(-2,2,200)
ax1.plot(strain1,Gc1,'ro',
         x,at.amplModel(par1,x),'r-')

ax1.legend(handles=[red_patch,g1_line,g2_dots,krit_line],loc='lower right',fontsize='14')
ax1.set_title('Amplitudni test - gelan',fontsize='18')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel("|G*| (Pa)",fontsize='16')
ax1.set_xlabel('$\gamma$ (%)', fontsize='16')
ax1.set_xlim((1.e-2,1.e2))
ax1.set_ylim((1,100))
ax1.axvline(krit1, color='r', linestyle='--')

ax12 = ax1.twinx()
ax12.plot(strain1,shear1,'rd')
ax12.set_yscale('log')
ax12.set_ylabel('$\\tau$ (Pa)', fontsize='16')
ax12.set_xlim((1.e-2,1.e2))

ax2 = plt.subplot(122)

x = np.logspace(-2,2,200)
ax2.plot(strain4,Gc4,'bo',
         x,at.amplModel(par4,x),'b-',
         strain5,Gc5,'ro',
         x,at.amplModel(par5,x),'r-')

ax2.legend(handles=[red_patch2,blue_patch2,g1_line,g2_dots,krit_line],loc='lower right',fontsize='12')
ax2.set_title('Amplitudni test - ksantan',fontsize='18')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel("|G*| (Pa)",fontsize='16')
ax2.set_xlabel('$\gamma$ (%)', fontsize='16')
ax2.set_xlim((1.e-2,1.e2))
ax2.set_ylim((1,100))
ax2.axvline(krit4, color='b', linestyle='--')
ax2.axvline(krit5, color='r', linestyle='--')

ax22 = ax2.twinx()
ax22.plot(strain4,shear4,'bd',
          strain5,shear5,'rd')
ax22.set_yscale('log')
ax22.set_ylabel('$\\tau$ (Pa)',fontsize='16')
ax22.set_xlim((1.e-2,1.e2))

plt.savefig("S2ampl.png")

#~ plt.show()
plt.close()


