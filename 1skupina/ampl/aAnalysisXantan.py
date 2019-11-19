import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import amplTest as at

# stevilo praznih vrstic
prazneVrstice = 5

# imena datotek s podatki
file1 = 'aS1XantanP1T20.csv'
file2 = 'aS1XantanP1T30.csv'

# uvozimo podatke
strain1, shear1, G1_1, G2_1 = at.podatki(file1, prazneVrstice)
strain2, shear2, G1_2, G2_2 = at.podatki(file2, prazneVrstice)

# izracunamo kompleksni modul
Gc1 = np.sqrt(G1_1**2+G2_1**2)
Gc2 = np.sqrt(G1_2**2+G2_2**2)

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

# natisnemo parametre na zaslon
print "k1", ["%0.9f" % i for i in k1]
print "k2", ["%0.9f" % i for i in k2]


#~ np.savetxt('S1_gelan1.txt',np.column_stack((g1,l1)),header='gelan 0.5 - g, l',delimiter=' ')
#~ np.savetxt('S1_gelan2.txt',np.column_stack((g2,l2)),header='gelan 0.7 - g, l',delimiter=' ')

# GRAF
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

blue_patch = mpatches.Patch(color='blue', label='Xantan 1.0 %, 20 $^\circ$C')
red_patch = mpatches.Patch(color='red', label='Xantan 1.0 %, 30 $^\circ$C')
g1_line = mlines.Line2D([], [], color='grey', marker='o',label="G*")

plt.figure(figsize=(15,5))
plt.gcf().subplots_adjust(bottom=0.15)
plt.subplots_adjust(hspace=1.0)

# MEHANSKI SPEKTER
plt.subplot(121)
x = np.logspace(-3,3,200)
plt.plot(strain1,Gc1,'bo',
         x,at.amplModel(par1,x),'b-',
         strain2,Gc2,'ro',
         x,at.amplModel(par2,x),'r-')

plt.legend(handles=[blue_patch,red_patch,g1_line],loc='lower left',fontsize='10')
plt.title('Amplitudni test',fontsize='18')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("G* (Pa)",fontsize='16')
plt.xlabel('$\gamma$ (%)', fontsize='16')

#~ # RELAKSACIJSKI SPEKTER
#~ plt.subplot(122)
#~ plt.plot(l1,g1,'ro-',
         #~ l2,g2,'bo-',
         #~ l3,g3,'go-')
#~ plt.legend(('gelan 0.5%', 'gelan 0.7%', 'gelan 1.0%'))
#~ plt.title('Relaksacijski spekter',fontsize='18')
#~ plt.yscale('log')
#~ plt.xscale('log')
#~ plt.ylabel("$g_i$ (Pa)",fontsize='16')
#~ plt.xlabel('$\lambda_i$ (s)',fontsize='16')
#~
#~ plt.savefig("S1_gelan.png")
plt.show()
#~ plt.close()


