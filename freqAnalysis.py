import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import maxwell as m
import freqTest as ft

# stevilo praznih vrstic
prazneVrstice = 5

# imena datotek s podatki
file1 = '1skupina.csv'

# uvozimo podatke
omega, G1, G2 = ft.podatki(file1, prazneVrstice)

# stevilo maxwellovih elementov
N = 5

# nastavimo zacetne vrednosti in omejitve parametrov
par = Parameters()
for i in range(1,N+1):
    gstr = 'g'+str(i)
    lstr = 'l'+str(i)
    par.add(gstr,value=1,min=1.e-4)
    par.add(lstr,value=1,min=1.e-4)

# minimizacija napake - non-linear least squares method
out = minimize(ft.maxwell,par,args=(N,omega,G1,G2),method='leastsq')

# pridobimo parametre iz fitanja
g, l = ft.koeficienti(par,N)

# natisnemo parametre na zaslon
print "g", ["%0.4f" % i for i in g]
print "l", ["%0.4f" % i for i in l]



# GRAF

# MEHANSKI SPEKTER
plt.subplot(121)
x = np.logspace(-3,8,200)
plt.plot(omega,G1,'o',
         omega,G2,'o',
         x,ft.maxwell(par,N,x)[0],
         x,ft.maxwell(par,N,x)[1])

plt.legend(('storage','loss','python S5','python L5'),loc='lower right')
plt.title('Mehanski spekter',fontsize='18')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("G', G'' (Pa)",fontsize='16')
plt.xlabel('$\omega$ (rad/s)', fontsize='16')

# RELAKSACIJSKI SPEKTER
plt.subplot(122)
plt.plot(l,g,'o')
plt.legend(('python'))
plt.title('Relaksacijski spekter',fontsize='18')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("$g_i$ (Pa)",fontsize='16')
plt.xlabel('$\lambda_i$ (s)',fontsize='16')

plt.show()

