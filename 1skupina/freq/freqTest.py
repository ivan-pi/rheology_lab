# Frekvencni test
# Definicije potrebne za analizo parametrov

import csv
import numpy as np
from lmfit import minimize, Parameters, fit_report

def podatki(filename,stVrst):

    # INITIALIZE LISTS
    omega = []
    storG = []
    lossG = []
    #~ dFactor = []
    #~ cVisc = []

    # READ FROM CSV FILE
    with open(filename, 'rb') as csvfile:
        for i in range(0,stVrst):
            next(csvfile)
        podatki = csv.reader(csvfile)
        for row in podatki:
            omega.append(float(row[1]))
            storG.append(float(row[2]))
            lossG.append(float(row[3]))
            #~ dFactor.append(float(row[4]))
            #~ cVisc.append(float(row[5]))

    # CONVERT LIST TO ARRAYS
    omega = np.array(omega)
    storG = np.array(storG)
    lossG = np.array(lossG)
    #~ dFactor = np.array(dFactor)
    #~ cVisc = np.array(cVisc)

    return omega, storG, lossG

def maxwell(params,N,omega,y1=None,y2=None):
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

def koeficienti(par,N):

    # EXTRACT PARAMETER VALUES
    p = par.valuesdict()
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

    return g, l;
