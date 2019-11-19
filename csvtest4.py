import csv
import numpy as np
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
import maxwell as m

# INITIALIZE LISTS
omega = []
storG = []
lossG = []
dFactor = []
cVisc = []

# READ FROM CSV FILE
with open('1skupina.csv', 'rb') as csvfile:
    for i in range(0,5):
        next(csvfile)
    podatki = csv.reader(csvfile)
    for row in podatki:
        omega.append(float(row[1]))
        storG.append(float(row[2]))
        lossG.append(float(row[3]))
        dFactor.append(float(row[4]))
        cVisc.append(float(row[5]))

# CONVERT LIST TO ARRAYS
omega = np.array(omega)
storG = np.array(storG)
lossG = np.array(lossG)
dFactor = np.array(dFactor)
cVisc = np.array(cVisc)

#~ p_true = Parameters()
#~ p_true.add('g1',value=10,min=0.0)
#~ p_true.add('g2',value=10,min=0.0)
#~ p_true.add('g3',value=1,min=0.0)
#~ p_true.add('g4',value=1,min=0.0)
#~ p_true.add('g5',value=1,min=0.0)
#~ p_true.add('l1',value=0.1,min=0.0)
#~ p_true.add('l2',value=1,min=0.0)
#~ p_true.add('l3',value=10,min=0.0)
#~ p_true.add('l4',value=100,min=0.0)
#~ p_true.add('l5',value=1000,min=0.0)

def residual(params,omega,y1=None,y2=None):
    g1 = params['g1'].value
    g2 = params['g2'].value
    g3 = params['g3'].value
    g4 = params['g4'].value
    g5 = params['g5'].value
    l1 = params['l1'].value
    l2 = params['l2'].value
    l3 = params['l3'].value
    l4 = params['l4'].value
    l5 = params['l5'].value

    sG = (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
         g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
         g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
         g4*l4**2*omega**2/(1.0+l4**2*omega**2)+
         g5*l5**2*omega**2/(1.0+l5**2*omega**2))

    lG = (g1*l1*omega/(1.0+l1**2*omega**2)+
         g2*l2*omega/(1.0+l2**2*omega**2)+
         g3*l3*omega/(1.0+l3**2*omega**2)+
         g4*l4*omega/(1.0+l4**2*omega**2)+
         g5*l5*omega/(1.0+l5**2*omega**2))

    if (y1 is None) and (y2 is None):
        return sG,lG

    return np.r_[y1-sG,
                 y2-lG];

params = Parameters()
params.add('g1',value=1,min=1.e-4)
params.add('g2',value=1,min=1.e-4)
params.add('g3',value=1,min=1.e-4)
params.add('g4',value=1,min=1.e-4)
params.add('g5',value=1,min=1.e-4)
params.add('l1',value=1,min=1.e-4)
params.add('l2',value=1,min=1.e-4)
params.add('l3',value=1,min=1.e-4)
params.add('l4',value=1,min=1.e-4)
params.add('l5',value=1,min=1.e-4)

out = minimize(residual,params,args=(omega,storG,lossG),method='leastsq')

fit = residual(params,omega,storG,lossG)
print(fit_report(params))

print fit

p = params.valuesdict()
g = np.array([p['g1'],p['g2'],p['g3'],p['g4'],p['g5']])
l = np.array([p['l1'],p['l2'],p['l3'],p['l4'],p['l5']])

plt.subplot(121)
x = np.logspace(-3,4,200)
plt.plot(omega,storG,'o',
         omega,lossG,'o',
         x,residual(params,x)[0],
         x,residual(params,x)[1])

plt.legend(('storage','loss','MS5','ML5'),loc='lower right')
plt.title('Frequency test')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("G',G''")
plt.xlabel('omega')

plt.subplot(122)
plt.plot(l,g,'o',)
plt.legend(('M5'))
plt.title('Relaxation specter')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("g_i")
plt.xlabel('lambda_i')

plt.show()

