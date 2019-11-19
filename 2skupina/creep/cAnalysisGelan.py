import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import creepTest as ct

# stevilo praznih vrstic
prazneVrstice = 5

# imena datotek s podatki
#~ file1 =  'cS1gelanP0.7T20.csv'
#~ file2 =  'cS1gelanP1.0T20.csv'
file3 = 'cS2xantanP1T20.csv'

# uvozimo podatke
#~ time1c,time1r, strain1c,strain1r, comp1 = ct.podatki(file1, prazneVrstice)
#~ time2c,time2r, strain2c,strain2r, comp2 = ct.podatki(file2, prazneVrstice)
time3c,time3r, strain3c,strain3r, comp3 = ct.podatki(file3, prazneVrstice)

comp3c = strain3c/0.1
comp3r = strain3r/0.1

# nastavimo zacetne vrednosti in omejitve parametrov
#~ par1 = Parameters()
#~ par1.add('G0' ,value=1.0,min=1.e-5,max=100)
#~ par1.add('l0',value=1.0,min=1.e-5,max=100)
#~ par1.add('G1' ,value=1.0,min=1.e-5,max=100)
#~ par1.add('l1' ,value=1.0,min=1.e-5,max=100)
#~ par1.add('tauC',value=0.03,vary=False)
#~
#~ par1kv = Parameters()
#~ par1kv.add('G' ,value=1.0,min=1.e-5)
#~ par1kv.add('l',value=1.0,min=1.e-5)
#~ par1kv.add('tauC',value=0.03,vary=False)
#~
#~ # minimizacija napake - non-linear least squares method
#~ out1 = minimize(ct.burgerModel,par1,args=(time1c,time1r,strain1c,strain1r),method='leastsq')
#~ out1kv = minimize(ct.KVModel,par1kv,args=(time1c,time1r,strain1c,strain1r),method='leastsq')
#~
#~ # pridobimo parametre iz fitanja
#~ k1 = ct.koefBurger(par1)
#~ k1kv = ct.koefKV(par1kv)
#~
#~ # nastavimo zacetne vrednosti in omejitve parametrov
#~ par2 = Parameters()
#~ par2.add('G0' ,value=1.0,min=1.e-5)
#~ par2.add('l0',value=1.0,min=1.e-5)
#~ par2.add('G1' ,value=1.0,min=1.e-5)
#~ par2.add('l1' ,value=1.0,min=1.e-5)
#~ par2.add('tauC',value=0.3,vary=False)
#~
#~ par2kv = Parameters()
#~ par2kv.add('G' ,value=1.0,min=1.e-5)
#~ par2kv.add('l',value=1.0,min=1.e-5)
#~ par2kv.add('tauC',value=0.3,vary=False)
#~
#~ # minimizacija napake - non-linear least squares method
#~ out2 = minimize(ct.burgerModel,par2,args=(time2c,time2r,strain2c,strain2r),method='leastsq')
#~ out2 = minimize(ct.KVModel,par2kv,args=(time2c,time2r,strain2c,strain2r),method='leastsq')
#~
#~ # pridobimo parametre iz fitanja
#~ k2 = ct.koefBurger(par2)
#~ k2kv = ct.koefKV(par2kv)

# nastavimo zacetne vrednosti in omejitve parametrov
par3 = Parameters()
par3.add('G0' ,value=10.0,min=1.e-5)
par3.add('l0',value=44.0,min=1.e-5)
par3.add('G1' ,value=11.0,min=1.e-5)
par3.add('l1' ,value=9.0,min=1.e-5)
par3.add('tauC',value=0.1,vary=False)

par3comp = Parameters()
par3comp.add('G0' ,value=10.0,min=1.e-5)
par3comp.add('l0',value=44.0,min=1.e-5)
par3comp.add('G1' ,value=11.0,min=1.e-5)
par3comp.add('l1' ,value=9.0,min=1.e-5)
par3comp.add('tauC',value=0.1,vary=False)

par3kv = Parameters()
par3kv.add('G' ,value=1.0,min=1.e-5)
par3kv.add('l',value=1.0,min=1.e-5)
par3kv.add('tauC',value=0.1,vary=False)

par3j = Parameters()
par3j.add('ni0' ,value=1.0,min=1.e-5)
par3j.add('G1' ,value=1.0,min=1.e-5)
par3j.add('l1' ,value=1.0,min=1.e-5)
par3j.add('tauC',value=0.1,vary=False)

# minimizacija napake - non-linear least squares method
out3 = minimize(ct.burgerModel,par3,args=(time3c,time3r,strain3c,strain3r),method='leastsq')
out3comp = minimize(ct.burgerModelJ,par3comp,args=(time3c,time3r,comp3c,comp3r),method='leastsq')
out3kv = minimize(ct.KVModel,par3kv,args=(time3c,time3r,strain3c,strain3r),method='leastsq')
out3j = minimize(ct.JeffreyModel,par3j,args=(time3c,time3r,strain3c,strain3r),method='leastsq')


# pridobimo parametre iz fitanja
k3 = ct.koefBurger(par3)
k3comp = ct.koefBurger(par3comp)
k3kv = ct.koefKV(par3kv)
k3j = ct.koefJeffrey(par3j)

# natisnemo parametre na zaslon
#~ print "k1","gelan 0.7", ["%0.12f" % i for i in k1]
#~ print "k2","gelan 1.0", ["%0.12f" % i for i in k2]
print "k3","xantan 1.0", ["%0.12f" % i for i in k3]

#~ print "k1","gelan 0.7", ["%0.12f" % i for i in k1]
#~ print "k2","gelan 1.0", ["%0.12f" % i for i in k2]
print "k3comp","xantan 1.0", ["%0.12f" % i for i in k3comp]

#~ print "k1kv","gelan 0.7", ["%0.12f" % i for i in k1kv]
#~ print "k2kv","gelan 1.0", ["%0.12f" % i for i in k2kv]
print "k3kv","xantan 1.0", ["%0.12f" % i for i in k3kv]

#~ print "k1j","gelan 0.7", ["%0.12f" % i for i in k1j]
#~ print "k2j","gelan 1.0", ["%0.12f" % i for i in k2j]
print "k3j","xantan 1.0", ["%0.12f" % i for i in k3j]

# GRAF
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import gridspec

red_patch = mpatches.Patch(color='red', label='gelan 0.7 %, 20 $^\circ$C')
blue_patch = mpatches.Patch(color='blue', label='gelan 1.0 %, 20 $^\circ$C')
green_patch = mlines.Line2D([], [], color='green',linestyle=' ',marker='o', label='xantan 1.0 %, 20 $^\circ$C')

g1_line = mlines.Line2D([], [], color='grey',linestyle=' ',marker='o',label="$\\gamma$")
g2_line = mlines.Line2D([], [], color='g',linestyle='-',label="Burger")
g3_line = mlines.Line2D([], [], color='g',linestyle='--',label="Jeffrey")
J_dots = mlines.Line2D([], [], color='grey', linestyle=' ',marker='o',label="J")
#~ krit_line = mlines.Line2D([], [], color='grey', linestyle='--',label="$\\gamma_{\mathrm{krit}}$")

plt.figure(figsize=(15,5))
plt.gcf().subplots_adjust(bottom=0.15)
plt.subplots_adjust(wspace=0.2)

# TEST LEZENJA IN OBNOVE
ax1 = plt.subplot(121)

xc = np.linspace(0,161,161)
xr = np.linspace(161,500,339)
ax1.plot(time3c,strain3c,'go',
         time3r,strain3r,'go',
         np.r_[xc,xr],ct.burgerModel(par3,xc,xr),'g-',
         np.r_[xc,xr],ct.JeffreyModel(par3j,xc,xr),'g--')
         #xc,ct.KVModel(par3kv,xc,xr)[0],'g--',
         #xr,ct.KVModel(par3kv,xc,xr)[1],'g--')

ax1.legend(handles=[green_patch,g2_line,g3_line],loc='lower right',fontsize='13',ncol=1)
ax1.set_title('Deformacija',fontsize='18')
ax1.set_ylabel('$\\gamma$ (/)',fontsize='16')
ax1.set_xlabel('$t$ (s)', fontsize='16')
ax1.set_xlim((0,350))
ax1.axvline(161, color='black', linestyle=':')


# VOLJNOST
ax2 = plt.subplot(122)

xc = np.linspace(0,161,161)
xr = np.linspace(161,500,339)
ax2.plot(time3c,comp3c,'go',
         time3r,comp3r,'go',
         np.r_[xc,xr],ct.burgerModelJ(par3comp,xc,xr),'g-')

ax2.legend(handles=[green_patch,g2_line],loc='lower right',fontsize='13')
ax2.set_title('Voljnost',fontsize='18')
ax2.set_ylabel("$J$ (Pa$^{-1}$)",fontsize='16')
ax2.set_xlabel('$t$ (s)', fontsize='16')
ax2.set_xlim((0,350))
ax2.axvline(161, color='black', linestyle=':')
plt.savefig("S2creep.png")
plt.show()
plt.close()


