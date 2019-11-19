import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import amplTest as at

# stevilo praznih vrstic
prazneVrstice = 5

# imena datotek s podatki
file1 = 'aS1GelanP0.5T20.csv'
file2 = 'aS1GelanP0.7T20.csv'
file3 = 'aS1GelanP1.0T20.csv'
file4 = 'aS1XantanP1T20.csv'
file5 = 'aS1XantanP1T30.csv'

# uvozimo podatke
strain1, shear1, G1_1, G2_1 = at.podatki(file1, prazneVrstice)
strain2, shear2, G1_2, G2_2 = at.podatki(file2, prazneVrstice)
strain3, shear3, G1_3, G2_3 = at.podatki(file3, prazneVrstice)
strain4, shear4, G1_4, G2_4 = at.podatki(file4, prazneVrstice)
strain5, shear5, G1_5, G2_5 = at.podatki(file5, prazneVrstice)

# izracunamo kompleksni modul
Gc1 = np.sqrt(G1_1**2+G2_1**2)
Gc2 = np.sqrt(G1_2**2+G2_2**2)
Gc3 = np.sqrt(G1_3**2+G2_3**2)
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
par2 = Parameters()
par2.add('G0',value=1,min=0.0)
par2.add('a' ,value=1,min=0.0)
par2.add('b' ,value=1,min=0.0)
par2.add('n' ,value=1,min=0.0)

# minimizacija napake - non-linear least squares method
out2 = minimize(at.amplModel,par2,args=(strain2,Gc2),method='leastsq')

# pridobimo parametre iz fitanja
k2 = at.koeficienti(par2)

# nastavimo zacetne vrednosti in omejitve parametrov
par3 = Parameters()
par3.add('G0',value=1,min=0.0)
par3.add('a' ,value=1,min=0.0)
par3.add('b' ,value=1,min=0.0)
par3.add('n' ,value=1,min=0.0)

# minimizacija napake - non-linear least squares method
out3 = minimize(at.amplModel,par3,args=(strain3,Gc3),method='leastsq')

# pridobimo parametre iz fitanja
k3 = at.koeficienti(par3)

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
par5.add('n' ,value=1,min=0.0)

# minimizacija napake - non-linear least squares method
out5 = minimize(at.amplModel,par5,args=(strain5,Gc5),method='leastsq')

# pridobimo parametre iz fitanja
k5 = at.koeficienti(par5)

krit1 = (0.03/(0.97*k1[1]-k1[2]))**(1.0/k1[3])
krit2 = (0.03/(0.97*k2[1]-k2[2]))**(1.0/k2[3])
krit3 = (0.03/(0.97*k3[1]-k3[2]))**(1.0/k3[3])
krit4 = (0.03/(0.97*k4[1]-k4[2]))**(1.0/k4[3])
krit5 = (0.03/(0.97*k5[1]-k5[2]))**(1.0/k5[3])

# natisnemo parametre na zaslon
print "kr1",krit1,"gelan 0.5", ["%0.12f" % i for i in k1]
print "kr2",krit2,"gelan 0.7", ["%0.12f" % i for i in k2]
print "kr3",krit3,"gelan 1.0", ["%0.12f" % i for i in k3]
print "kr4",krit4,"xantan 20C", ["%0.12f" % i for i in k4]
print "kr5",krit5,"xantan 30C", ["%0.12f" % i for i in k5]

#~ np.savetxt('S1_gelan1.txt',np.column_stack((g1,l1)),header='gelan 0.5 - g, l',delimiter=' ')
#~ np.savetxt('S1_gelan2.txt',np.column_stack((g2,l2)),header='gelan 0.7 - g, l',delimiter=' ')
#~ np.savetxt('S1_gelan3.txt',np.column_stack((g3,l3)),header='gelan 1.0 - g, l',delimiter=' ')

# GRAF
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

red_patch = mpatches.Patch(color='red', label='gelan 0.5 %, 20 $^\circ$C')
blue_patch = mpatches.Patch(color='blue', label='gelan 0.7 %, 20 $^\circ$C')
green_patch = mpatches.Patch(color='green', label='gelan 1.0 %, 20 $^\circ$C')
blue_patch2 = mpatches.Patch(color='blue', label='ksantan 1 %, 20 $^\circ$C')
red_patch2 = mpatches.Patch(color='red', label='ksantan 1 %, 30 $^\circ$C')
g1_line = mlines.Line2D([], [], color='grey', marker='o',label="G*")
tau_dots = mlines.Line2D([], [], color='grey', marker='d',label="$\\tau$")
krit_line = mlines.Line2D([], [], color='grey', linestyle='--',label="$\\gamma_{\mathrm{krit}}$")

plt.figure(figsize=(15,5))
plt.gcf().subplots_adjust(bottom=0.15)
plt.subplots_adjust(wspace=0.4)

# MEHANSKI SPEKTER
ax1 = plt.subplot(121)
x = np.logspace(-2,2,200)

ax1.plot(strain1,Gc1,'ro',
         x,at.amplModel(par1,x),'r-',
         strain2,Gc2,'bo',
         x,at.amplModel(par2,x),'b-',
         strain3,Gc3,'go',
         x,at.amplModel(par3,x),'g-')

ax1.legend(handles=[red_patch,blue_patch,green_patch,g1_line,tau_dots,krit_line],loc='lower right',fontsize='10',ncol=2)
ax1.set_title('Amplitudni test - gelan',fontsize='18')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel("|G*| (Pa)",fontsize='16')
ax1.set_xlabel('$\gamma$ (%)', fontsize='16')
ax1.axvline(krit1, color='r', linestyle='--')
ax1.axvline(krit2, color='b', linestyle='--')
ax1.axvline(krit3, color='g', linestyle='--')
ax1.set_xlim((0.01,100))
ax1.set_ylim((1,100))

ax12 = ax1.twinx()
ax12.plot(strain1,shear1,'rd',
          strain2,shear2,'bd',
          strain3,shear3,'gd')
ax12.set_yscale('log')
ax12.set_ylabel('$\\tau$ (Pa)', fontsize='16')
ax12.set_xlim((0.01,100))

ax2 = plt.subplot(122)

x = np.logspace(-2,2,200)
ax2.plot(strain4,Gc4,'bo',
         x,at.amplModel(par4,x),'b-',
         strain5,Gc5,'ro',
         x,at.amplModel(par5,x),'r-')

ax2.legend(handles=[blue_patch2,red_patch2,g1_line,tau_dots,krit_line],loc='lower right',fontsize='12')
ax2.set_title('Amplitudni test - ksantan',fontsize='18')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel("|G*| (Pa)",fontsize='16')
ax2.set_xlabel('$\gamma$ (%)', fontsize='16')
ax2.axvline(krit4, color='b', linestyle='--')
ax2.axvline(krit5, color='r', linestyle='--')
ax2.set_xlim((0.01,100))
ax2.set_ylim((1,100))

ax22 = ax2.twinx()
ax22.plot(strain4,shear4,'bd',
          strain5,shear5,'rd')
ax22.set_yscale('log')
ax22.set_ylabel('$\\tau$ (Pa)', fontsize='16')
ax22.set_xlim((0.01,100))
ax22.set_ylim((0.0001,100))

plt.savefig("S1ampl.png")
#~ plt.show()
plt.close()


