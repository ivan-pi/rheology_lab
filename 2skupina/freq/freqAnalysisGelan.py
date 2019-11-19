import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import maxwell as m
import freqTest as ft

# stevilo praznih vrstic
prazneVrstice = 5

# imena datotek s podatki
file1 = 'fS2gelanP0.5T20.csv'
file2 = 'fS2gelanP0.7T20.csv'
file3 = 'fS2gelanP1.0T20.csv'

# uvozimo podatke
omega1, G1_1, G2_1 = ft.podatki(file1, prazneVrstice)
omega2, G1_2, G2_2 = ft.podatki(file2, prazneVrstice)
omega3, G1_3, G2_3 = ft.podatki(file3, prazneVrstice)

# stevilo maxwellovih elementov
N1 = 6
N2 = 6
N3 = 6

# nastavimo zacetne vrednosti in omejitve parametrov
par1 = Parameters()
for i in range(1,N1+1):
    gstr = 'g'+str(i)
    lstr = 'l'+str(i)
    par1.add(gstr,value=1,min=1.e-2)
    par1.add(lstr,value=1,min=1.e-4)

# minimizacija napake - non-linear least squares method
out1 = minimize(ft.maxwell,par1,args=(N1,omega1,G1_1,G2_1),method='leastsq')

# pridobimo parametre iz fitanja
g1, l1 = ft.koeficienti(par1,N1)

# nastavimo zacetne vrednosti in omejitve parametrov
par2 = Parameters()
for i in range(1,N2+1):
    gstr = 'g'+str(i)
    lstr = 'l'+str(i)
    par2.add(gstr,value=1,min=1.e-4)
    par2.add(lstr,value=1,min=1.e-4)

# minimizacija napake - non-linear least squares method
out2 = minimize(ft.maxwell,par2,args=(N2,omega2,G1_2,G2_2),method='leastsq')

# pridobimo parametre iz fitanja
g2, l2 = ft.koeficienti(par2,N2)

# nastavimo zacetne vrednosti in omejitve parametrov
par3 = Parameters()
for i in range(1,N3+1):
    gstr = 'g'+str(i)
    lstr = 'l'+str(i)
    par3.add(gstr,value=1,min=1.e-4)
    par3.add(lstr,value=1,min=1.e-4)

# minimizacija napake - non-linear least squares method
out3 = minimize(ft.maxwell,par3,args=(N3,omega3,G1_3,G2_3),method='leastsq')

# pridobimo parametre iz fitanja
g3, l3 = ft.koeficienti(par3,N3)

# natisnemo parametre na zaslon
print "g1", ["%0.4f" % i for i in g1]
print "l1", ["%0.4f" % i for i in l1]

# natisnemo parametre na zaslon
print "g2", ["%0.4f" % i for i in g2]
print "l2", ["%0.4f" % i for i in l2]

# natisnemo parametre na zaslon
print "g3", ["%0.4f" % i for i in g3]
print "l3", ["%0.4f" % i for i in l3]

#~ np.savetxt('S1_gelan1.txt',np.column_stack((g1,l1)),header='gelan 0.5 - g, l',delimiter=' ')
#~ np.savetxt('S1_gelan2.txt',np.column_stack((g2,l2)),header='gelan 0.7 - g, l',delimiter=' ')
#~ np.savetxt('S1_gelan3.txt',np.column_stack((g3,l3)),header='gelan 1.0 - g, l',delimiter=' ')

# GRAF
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

red_patch = mpatches.Patch(color='red', label='gelan 0.5 %')
blue_patch = mpatches.Patch(color='blue', label='gelan 0.7 %')
green_patch = mpatches.Patch(color='green', label='gelan 1.0%')
g1_line = mlines.Line2D([], [], color='grey', linestyle='-',  marker='o',label="G'")
g2_line = mlines.Line2D([], [], color='grey', linestyle='--', marker='d',label="G''")

plt.figure(figsize=(15,5))
plt.gcf().subplots_adjust(bottom=0.15)
plt.subplots_adjust(hspace=1.0)

# MEHANSKI SPEKTER
plt.subplot(121)
x = np.logspace(-2,3,200)
plt.plot(omega1,G1_1,'ro',
         omega1,G2_1,'rd',
         x,ft.maxwell(par1,N1,x)[0],'r-',
         x,ft.maxwell(par1,N1,x)[1],'r--',
         omega2,G1_2,'bo',
         omega2,G2_2,'bd',
         x,ft.maxwell(par2,N2,x)[0],'b-',
         x,ft.maxwell(par2,N2,x)[1],'b--',
         omega3,G1_3,'go',
         omega3,G2_3,'gd',
         x,ft.maxwell(par3,N3,x)[0],'g-',
         x,ft.maxwell(par3,N3,x)[1],'g--')

plt.legend(("gelan 0.5% G'",
            "gelan 0.5% G''",
            "model G'",
            "model G''"),loc='lower right')

plt.legend(handles=[red_patch,blue_patch,green_patch,g1_line,g2_line],loc='lower right')
plt.title('Mehanski spekter',fontsize='18')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("G', G'' (Pa)",fontsize='16')
plt.xlabel('$\omega$ (rad/s)', fontsize='16')

# RELAKSACIJSKI SPEKTER
plt.subplot(122)
plt.plot(l1,g1,'ro-',
         l2,g2,'bo-',
         l3,g3,'go-')
plt.legend(('gelan 0.5%', 'gelan 0.7%', 'gelan 1.0%'))
plt.title('Relaksacijski spekter',fontsize='18')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("$g_i$ (Pa)",fontsize='16')
plt.xlabel('$\lambda_i$ (s)',fontsize='16')

plt.savefig("S2_gelan.png")
#plt.show()
plt.close()


