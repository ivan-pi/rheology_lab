import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import creepTest as ct

# stevilo praznih vrstic
prazneVrstice = 5

# imena datotek s podatki
file1 =  'cS1gelanP0.7T20.csv'
file2 =  'cS1gelanP1.0T20.csv'
file3 = 'cS1xantanP1.0T20.csv'

# uvozimo podatke
time1c,time1r, strain1c,strain1r, comp1 = ct.podatki(file1, prazneVrstice)
time2c,time2r, strain2c,strain2r, comp2 = ct.podatki(file2, prazneVrstice)
time3c,time3r, strain3c,strain3r, comp3 = ct.podatki(file3, prazneVrstice)

comp1c = strain1c/0.03
comp1r = strain1r/0.03
comp2c = strain2c/0.3
comp2r = strain2r/0.3
comp3c = strain3c/0.1
comp3r = strain3r/0.1

# nastavimo zacetne vrednosti in omejitve parametrov
par1 = Parameters()
par1.add('G0' ,value=1.0,min=1.e-5,max=100)
par1.add('l0',value=1.0,min=1.e-5,max=100)
par1.add('G1' ,value=1.0,min=1.e-5,max=100)
par1.add('l1' ,value=1.0,min=1.e-5,max=100)
par1.add('tauC',value=0.03,vary=False)

par1comp = Parameters()
par1comp.add('G0' ,value=1.0,min=1.e-5,max=100)
par1comp.add('l0',value=1.0,min=1.e-5,max=100)
par1comp.add('G1' ,value=1.0,min=1.e-5,max=100)
par1comp.add('l1' ,value=1.0,min=1.e-5,max=100)
par1comp.add('tauC',value=0.03,vary=False)

par1j = Parameters()
par1j.add('ni0' ,value=100.0,min=1.e-5)
par1j.add('G1' ,value=0.224,min=1.e-5)
par1j.add('l1' ,value=45.373,min=1.e-5)
par1j.add('tauC',value=0.03,vary=False)

par1kv = Parameters()
par1kv.add('G' ,value=1.0,min=1.e-5)
par1kv.add('l',value=1.0,min=1.e-5)
par1kv.add('tauC',value=0.03,vary=False)

# minimizacija napake - non-linear least squares method
out1 = minimize(ct.burgerModel,par1,args=(time1c,time1r,strain1c,strain1r),method='leastsq')
out1comp = minimize(ct.burgerModelJ,par1comp,args=(time1c,time1r,comp1c,comp1r),method='leastsq')
out1j = minimize(ct.JeffreyModel,par1j,args=(time1c,time1r,strain1c,strain1r),method='leastsq')
out1kv = minimize(ct.KVModel,par1kv,args=(time1c,time1r,strain1c,strain1r),method='leastsq')

# pridobimo parametre iz fitanja
k1 = ct.koefBurger(par1)
k1comp = ct.koefBurger(par1comp)
k1j = ct.koefJeffrey(par1j)
k1kv = ct.koefKV(par1kv)

# nastavimo zacetne vrednosti in omejitve parametrov
par2 = Parameters()
par2.add('G0' ,value=1.0,min=1.e-5)
par2.add('l0',value=1.0,min=1.e-5)
par2.add('G1' ,value=1.0,min=1.e-5)
par2.add('l1' ,value=1.0,min=1.e-5)
par2.add('tauC',value=0.3,vary=False)

par2comp = Parameters()
par2comp.add('G0' ,value=1.0,min=1.e-5,max=100)
par2comp.add('l0',value=1.0,min=1.e-5,max=100)
par2comp.add('G1' ,value=1.0,min=1.e-5,max=100)
par2comp.add('l1' ,value=1.0,min=1.e-5,max=100)
par2comp.add('tauC',value=0.3,vary=False)

par2j = Parameters()
par2j.add('ni0' ,value=1.0,min=1.e-5)
par2j.add('G1' ,value=1.0,min=1.e-5)
par2j.add('l1' ,value=1.0,min=1.e-5)
par2j.add('tauC',value=0.3,vary=False)

par2kv = Parameters()
par2kv.add('G' ,value=1.0,min=1.e-5)
par2kv.add('l',value=1.0,min=1.e-5)
par2kv.add('tauC',value=0.3,vary=False)

# minimizacija napake - non-linear least squares method
out2 = minimize(ct.burgerModel,par2,args=(time2c,time2r,strain2c,strain2r),method='leastsq')
out2comp = minimize(ct.burgerModelJ,par2comp,args=(time2c,time2r,comp2c,comp2r),method='leastsq')
out2j = minimize(ct.JeffreyModel,par2j,args=(time2c,time2r,strain2c,strain2r),method='leastsq')
out2 = minimize(ct.KVModel,par2kv,args=(time2c,time2r,strain2c,strain2r),method='leastsq')

# pridobimo parametre iz fitanja
k2 = ct.koefBurger(par2)
k2comp = ct.koefBurger(par2comp)
k2j = ct.koefJeffrey(par2j)
k2kv = ct.koefKV(par2kv)

# nastavimo zacetne vrednosti in omejitve parametrov
par3 = Parameters()
par3.add('G0' ,value=1.0,min=1.e-5)
par3.add('l0',value=1.0,min=1.e-5)
par3.add('G1' ,value=1.0,min=1.e-5)
par3.add('l1' ,value=1.0,min=1.e-5)
par3.add('tauC',value=0.1,vary=False)

par3comp = Parameters()
par3comp.add('G0' ,value=1.0,min=1.e-5,max=100)
par3comp.add('l0',value=1.0,min=1.e-5,max=100)
par3comp.add('G1' ,value=1.0,min=1.e-5,max=100)
par3comp.add('l1' ,value=1.0,min=1.e-5,max=100)
par3comp.add('tauC',value=0.1,vary=False)

par3j = Parameters()
par3j.add('ni0' ,value=1.0,min=1.e-5)
par3j.add('G1' ,value=1.0,min=1.e-5)
par3j.add('l1' ,value=1.0,min=1.e-5)
par3j.add('tauC',value=0.1,vary=False)

par3kv = Parameters()
par3kv.add('G' ,value=1.0,min=1.e-5)
par3kv.add('l',value=1.0,min=1.e-5)
par3kv.add('tauC',value=0.1,vary=False)

# minimizacija napake - non-linear least squares method
out3 = minimize(ct.burgerModel,par3,args=(time3c,time3r,strain3c,strain3r),method='leastsq')
out3comp = minimize(ct.burgerModelJ,par3comp,args=(time3c,time3r,comp3c,comp3r),method='leastsq')
out3j = minimize(ct.JeffreyModel,par3j,args=(time3c,time3r,strain3c,strain3r),method='leastsq')
out3kv = minimize(ct.KVModel,par3kv,args=(time3c,time3r,strain3c,strain3r),method='leastsq')

# pridobimo parametre iz fitanja
k3 = ct.koefBurger(par3)
k3comp = ct.koefBurger(par3comp)
k3j = ct.koefJeffrey(par3j)
k3kv = ct.koefKV(par3kv)

# natisnemo parametre na zaslon
print "k1","gelan 0.7", ["%0.8f" % i for i in k1]
print "k2","gelan 1.0", ["%0.8f" % i for i in k2]
print "k3","xantan 1.0", ["%0.8f" % i for i in k3]

print "k1comp","gelan 0.7", ["%0.8f" % i for i in k1comp]
print "k2comp","gelan 1.0", ["%0.8f" % i for i in k2comp]
print "k3comp","xantan 1.0", ["%0.8f" % i for i in k3comp]

print "k1j","gelan 0.7", ["%0.8f" % i for i in k1j]
print "k2j","gelan 1.0", ["%0.8f" % i for i in k2j]
print "k3j","xantan 1.0", ["%0.8f" % i for i in k3j]

print "k1kv","gelan 0.7", ["%0.8f" % i for i in k1kv]
print "k2kv","gelan 1.0", ["%0.8f" % i for i in k2kv]
print "k3kv","xantan 1.0", ["%0.8f" % i for i in k3kv]
# GRAF
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import gridspec

red_patch = mpatches.Patch(color='red', label='gelan 0.7 %, 20 $^\circ$C')
blue_patch = mpatches.Patch(color='blue', label='gelan 1.0 %, 20 $^\circ$C')
green_patch = mpatches.Patch(color='green', label='xantan 1.0 %, 20 $^\circ$C')

g1_line = mlines.Line2D([], [], color='grey',linestyle=' ',marker='o',label="$\\gamma$")
g2_line = mlines.Line2D([], [], color='grey',linestyle='-',label="Burger")
g3_line = mlines.Line2D([], [], color='grey',linestyle='--',label="Jeffrey")
J_dots = mlines.Line2D([], [], color='grey', linestyle=' ',marker='o',label="J")
#~ krit_line = mlines.Line2D([], [], color='grey', linestyle='--',label="$\\gamma_{\mathrm{krit}}$")

plt.figure(figsize=(15,5))
plt.gcf().subplots_adjust(bottom=0.15)
plt.subplots_adjust(wspace=0.2)

# TEST LEZENJA IN OBNOVE
ax1 = plt.subplot(121)

xc = np.linspace(0,161,161)
xr = np.linspace(161,350,189)
ax1.plot(time1c,strain1c,'ro',
         time1r,strain1r,'ro',
         np.r_[xc,xr],ct.burgerModel(par1,xc,xr),'r-',
         np.r_[xc,xr],ct.JeffreyModel(par1j,xc,xr),'r--',
         time2c,strain2c,'bo',
         time2r,strain2r,'bo',
         np.r_[xc,xr],ct.burgerModel(par2,xc,xr),'b-',
         np.r_[xc,xr],ct.JeffreyModel(par2j,xc,xr),'b--',
         time3c,strain3c,'go',
         time3r,strain3r,'go',
         np.r_[xc,xr],ct.burgerModel(par3,xc,xr),'g-',
         np.r_[xc,xr],ct.JeffreyModel(par3j,xc,xr),'g--',)

ax1.legend(handles=[red_patch,blue_patch,green_patch,g1_line,g2_line,g3_line],loc='upper right',fontsize='11',ncol=1)
ax1.set_title('Deformacija',fontsize='18')
ax1.set_ylabel('$\\gamma$ (/)',fontsize='16')
ax1.set_xlabel('$t$ (s)', fontsize='16')
ax1.axvline(161, color='black', linestyle=':')


# VOLJNOST
ax2 = plt.subplot(122)

xc = np.linspace(0,161,161)
xr = np.linspace(161,350,189)
ax2.plot(time1c,comp1c,'ro',
         time1r,comp1r,'ro',
         np.r_[xc,xr],ct.burgerModelJ(par1comp,xc,xr),'r-',
         time2c,comp2c,'bo',
         time2r,comp2r,'bo',
         np.r_[xc,xr],ct.burgerModelJ(par2comp,xc,xr),'b-',
         time3c,comp3c,'go',
         time3r,comp3r,'go',
         np.r_[xc,xr],ct.burgerModelJ(par3comp,xc,xr),'g-',)

ax2.legend(handles=[red_patch,blue_patch,green_patch,J_dots,g2_line],loc='upper right',fontsize='11')
ax2.set_title('Voljnost',fontsize='18')
ax2.set_ylabel("$J$ (Pa$^{-1}$)",fontsize='16')
ax2.set_xlabel('$t$ (s)', fontsize='16')
ax2.axvline(161, color='black', linestyle=':')
plt.savefig("S1creep.png")
plt.show()
plt.close()


