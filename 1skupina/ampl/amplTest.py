# Amplitudni test
# Definicije potrebne za analizo parametrov

import csv
import numpy as np
from lmfit import minimize, Parameters, fit_report

def podatki(filename,stVrst):

    # INITIALIZE LISTS
    strain = []
    shearStress = []
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
            strain.append(float(row[1]))
            shearStress.append(float(row[2]))
            storG.append(float(row[3]))
            lossG.append(float(row[4]))
            #~ dFactor.append(float(row[4]))
            #~ cVisc.append(float(row[5]))

    # CONVERT LIST TO ARRAYS
    strain = np.array(strain)
    shearStress = np.array(shearStress)
    storG = np.array(storG)
    lossG = np.array(lossG)
    #~ dFactor = np.array(dFactor)
    #~ cVisc = np.array(cVisc)

    return strain, shearStress, storG, lossG

def amplModel(params,x,y=None):
    a = params['a'].value
    b = params['b'].value
    n = params['n'].value
    G0 = params['G0'].value

    Gcalc = G0*(1.0+b*x**n)/(1.0+a*x**n)

    if (y is None):
        return Gcalc

    return Gcalc-y;

def koeficienti(par):

    # EXTRACT PARAMETER VALUES
    p = par.valuesdict()
    koef = np.r_[p['G0'],p['a'],p['b'],p['n']]

    return koef;
