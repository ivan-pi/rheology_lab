import csv
import numpy as np
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


# Maxwell 1
stor1, scov1 = curve_fit(m.sG1fit,omega,storG)
print "stor1", stor1
loss1, lcov1 = curve_fit(m.lG1fit,omega,lossG)
print "loss1", loss1

# Maxwell 2
stor2, scov2 = curve_fit(m.sG2fit,omega,storG)
print "stor2", stor2
loss2, lcov2 = curve_fit(m.lG2fit,omega,lossG)
print "loss2", loss2

# Maxwell 3
#~ stor3, scov3 = curve_fit(m.sG3fit,omega,storG)
#~ print stor3
#~ loss3, lcov3 = curve_fit(m.lG3fit,omega,lossG)
#~ print loss3

# Maxwell 4
stor4, scov4 = curve_fit(m.sG4fit,omega,storG)
print "stor4", stor4
loss4, lcov4 = curve_fit(m.lG4fit,omega,lossG)
print "loss4", loss4

# Maxwell 5
stor5, scov5 = curve_fit(m.sG5fit,omega,storG)
print "stor5", stor5
loss5, lcov5 = curve_fit(m.lG5fit,omega,lossG)
print "loss5", loss5

# Maxwell 6
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
plt.plot(omega,storG,'o',
         omega,lossG,'o',
         x,m.sG1eval(x,stor1),x,m.lG1eval(x,loss1),
         x,m.sG2eval(x,stor2),x,m.lG2eval(x,loss2),
         #~ x,m.sG3eval(x,stor3),x,m.lG3eval(x,loss3),
         x,m.sG4eval(x,stor4),x,m.lG4eval(x,loss4),
         x,m.sG5eval(x,stor5),x,m.lG5eval(x,loss5))
         #~ x,m.sG6eval(x,stor6),x,m.lG6eval(x,loss6),
         #~ x,m.sG7eval(x,stor7),x,m.lG7eval(x,loss7),
         #~ x,m.sG8eval(x,stor8),x,m.lG8eval(x,loss8))

plt.title('Frequency test')
plt.yscale('log')
plt.xscale('log')
plt.ylabel("G',G''")
plt.xlabel('omega')
plt.show()




