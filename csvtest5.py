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

N = 5

def residual(params,N,omega,y1=None,y2=None):
    g = [None]*N
    l = [None]*N
    for i in range(1,N+1):
        gstr = 'g'+str(i)
        lstr = 'l'+str(i)
        g[i-1] = params[gstr].value
        l[i-1] = params[lstr].value

    sG = 0
    lG = 0
    for i in range(1,N+1):
        sG = sG + g[i-1]*l[i-1]**2*omega**2/(1.0+l[i-1]**2*omega**2)
        lG = lG + g[i-1]*l[i-1]*omega/(1.0+l[i-1]**2*omega**2)

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
#~ params.add('g6',value=1,min=1.e-4)
#~ params.add('g7',value=1,min=1.e-4)
#~ params.add('g8',value=1,min=1.e-4)
#~ params.add('g9',value=1,min=1.e-4)
#~ params.add('g10',value=1,min=1.e-4)
params.add('l1',value=1,min=1.e-4)
params.add('l2',value=1,min=1.e-4)
params.add('l3',value=1,min=1.e-4)
params.add('l4',value=1,min=1.e-4)
params.add('l5',value=1,min=1.e-4)
#~ params.add('l6',value=1,min=1.e-4)
#~ params.add('l7',value=1,min=1.e-4)
#~ params.add('l8',value=1,min=1.e-4)
#~ params.add('l9',value=1,min=1.e-4)
#~ params.add('l10',value=1,min=1.e-4)

out = minimize(residual,params,args=(N,omega,storG,lossG),method='leastsq')

fit = residual(params,N,omega)[0]
print "fit", fit
print(fit_report(params))

# EXTRACT PARAMETER VALUES
p = params.valuesdict()
g = np.array([])
l = np.array([])
for i in range(1,N+1):
    gstr = 'g'+str(i)
    lstr = 'l'+str(i)
    g = np.append(g,p[gstr])
    l = np.append(l,p[lstr])

# sortiramo v pravem zaporedju
g = g[np.argsort(l)]
l.sort()


print "g", np.round(g,4)
print "l", np.round(l,4)

g2 = np.r_[907.6350,10.6665,6.8126,5.6231,13.8406]
l2 = np.r_[0.0002,0.0387,0.2233,1.2070,17.8452]

# MEHANSKI SPEKTER
plt.subplot(121)
x = np.logspace(-3,4,200)
plt.plot(omega,storG,'o',
         omega,lossG,'o',
         x,residual(params,N,x)[0],
         x,residual(params,N,x)[1],
         x,m.sG5eval(x,np.append(g2,l2)),
         x,m.lG5eval(x,np.append(g2,l2)))

plt.legend(('storage','loss','python S5','python L5','MATLAB S5','MATLAB L5'),loc='lower right')
plt.title('Mehanski spekter',fontsize='18')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("G', G'' (Pa)",fontsize='16')
plt.xlabel('$\omega$ (rad/s)', fontsize='16')

# RELAKSACIJSKI SPEKTER
plt.subplot(122)
plt.plot(l,g,
         l2,g2,'o')
plt.legend(('python','MATLAB'))
plt.title('Relaksacijski spekter',fontsize='18')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("$g_i$ (Pa)",fontsize='16')
plt.xlabel('$\lambda_i$ (s)',fontsize='16')

plt.show()

