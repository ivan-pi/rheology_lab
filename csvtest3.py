import csv
import numpy as np
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import maxwell as m

omega = []
storG = []
lossG = []
dFactor = []
cVisc = []

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

omega = np.array(omega)
storG = np.array(storG)
lossG = np.array(lossG)
dFactor = np.array(dFactor)
cVisc = np.array(cVisc)
# 4 Maxwell elements STORAGE modulus
def sGfit4(omega,par):
    g1,g2,g3,g4,l1,l2,l3,l4 = par
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2));

# 4 Maxwell elements LOSS modulus
def lGfit4(omega,par):
    g1,g2,g3,g4,l1,l2,l3,l4 = par
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2));

def errFunc4(p,x,y1,y2):
    return np.r_[sGfit4(x,p)-y1,
                 lGfit4(x,p)-y2];

# 5 Maxwell elements STORAGE modulus
def sGfit5(omega,par):
    g1,g2,g3,g4,g5,l1,l2,l3,l4,l5 = par
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2)+
            g5*l5**2*omega**2/(1.0+l5**2*omega**2));

# 5 Maxwell elements LOSS modulus
def lGfit5(omega,par):
    g1,g2,g3,g4,g5,l1,l2,l3,l4,l5 = par
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2)+
            g5*l5*omega/(1.0+l5**2*omega**2));

def errFunc5(p,x,y1,y2):
    return np.r_[sGfit5(x,p)-y1,
                 lGfit5(x,p)-y2];

# 6 Maxwell elements STORAGE modulus
def sGfit6(omega,par):
    g1,g2,g3,g4,g5,g6,l1,l2,l3,l4,l5,l6 = par
    return (g1*l1**2*omega**2/(1.0+l1**2*omega**2)+
            g2*l2**2*omega**2/(1.0+l2**2*omega**2)+
            g3*l3**2*omega**2/(1.0+l3**2*omega**2)+
            g4*l4**2*omega**2/(1.0+l4**2*omega**2)+
            g5*l5**2*omega**2/(1.0+l5**2*omega**2)+
            g6*l6**2*omega**2/(1.0+l6**2*omega**2));

# 6 Maxwell elements LOSS modulus
def lGfit6(omega,par):
    g1,g2,g3,g4,g5,g6,l1,l2,l3,l4,l5,l6 = par
    return (g1*l1*omega/(1.0+l1**2*omega**2)+
            g2*l2*omega/(1.0+l2**2*omega**2)+
            g3*l3*omega/(1.0+l3**2*omega**2)+
            g4*l4*omega/(1.0+l4**2*omega**2)+
            g5*l5*omega/(1.0+l5**2*omega**2)+
            g6*l6*omega/(1.0+l6**2*omega**2));

def errFunc6(p,x,y1,y2):
    return np.r_[sGfit6(x,p)-y1,
                 lGfit6(x,p)-y2];

p6 = np.array([10,10,10,1.0,1.0,1.0,0.01,0.1,1,10.0,100,1000])
p_best6, pcov = leastsq(errFunc6,p6,args=(omega,storG,lossG),maxfev=20000)

p5 = np.array([10,10,10,1.0,1.0,0.01,0.1,1,10.0,100])
p_best5, pcov = leastsq(errFunc5,p5,args=(omega,storG,lossG),maxfev=20000)

p4 = np.array([10,10,1.0,1.0,0.1,1,10.0,100])
p_best4, pcov = leastsq(errFunc4,p4,args=(omega,storG,lossG))

print p_best6
print p_best5
print p_best4

#~ # Maxwell 3
#~ stor3, scov3 = curve_fit(m.sG3fit,omega,storG)
#~ print stor3
#~ loss3, lcov3 = curve_fit(m.lG3fit,omega,lossG)
#~ print loss3
#~
#~ # Maxwell 4
#~ stor4, scov4 = curve_fit(m.sG4fit,omega,storG)
#~ print "stor4", stor4
#~ loss4, lcov4 = curve_fit(m.lG4fit,omega,lossG)
#~ print "loss4", loss4
#~
#~ # Maxwell 5
#~ stor5, scov5 = curve_fit(m.sG5fit,omega,storG)
#~ print "stor5", stor5
#~ loss5, lcov5 = curve_fit(m.lG5fit,omega,lossG)
#~ print "loss5", loss5
#~
#~ # Maxwell 6
#~ stor6, scov6 = curve_fit(m.sG6fit,omega,storG)
#~ print stor6
#~ loss6, lcov6 = curve_fit(m.lG6fit,omega,lossG)
#~ print loss6

# Maxwell 7
#~ stor7, scov7 = curve_fit(m.sG7fit,omega,storG)
#~ print stor7
#~ loss7, lcov7 = curve_fit(m.lG7fit,omega,lossG)
#~ print loss7

# Maxwell 8
#~ stor8, scov8 = curve_fit(m.sG8fit,omega,storG)
#~ print stor8
#~ loss8, lcov8 = curve_fit(m.lG8fit,omega,lossG)
#~ print loss8

x = np.logspace(-1,3,100)
plt.figure(1)
plt.subplot(121)
plt.plot(omega,storG,'o',
         omega,lossG,'o',
         #~ x,m.sG4eval(x,stor4),
         #~ x,m.lG4eval(x,loss4),
         #~ x,m.sG5eval(x,stor5),
         #~ x,m.lG5eval(x,loss5),
         #~ x,m.sG6eval(x,stor6),
         #~ x,m.lG6eval(x,loss6),
         x,m.sG4eval(x,p_best4),
         x,m.lG4eval(x,p_best4),
         x,m.sG5eval(x,p_best5),
         x,m.lG5eval(x,p_best5),
         x,m.sG6eval(x,p_best6),
         x,m.lG6eval(x,p_best6))

plt.legend(('storage','loss','SS4','LL4','SS5','LL5','SS6','LL6'))
plt.title('Frequency test')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("G',G''")
plt.xlabel('omega')

plt.subplot(122)
plt.plot(p_best4[0:3],p_best4[4:7],'o',
         p_best5[0:4],p_best5[5:9],'o',
         p_best6[0:5],p_best6[6:11],'o')

plt.legend(('M4','M5','M6'))
plt.title('Relaxation specter')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("g_i")
plt.xlabel('lambda_i')

plt.show()

