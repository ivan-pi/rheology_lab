# Amplitudni test
# Definicije potrebne za analizo parametrov

import csv
import numpy as np
from lmfit import minimize, Parameters, fit_report

def podatki(filename,stVrst):

    # INITIALIZE LISTS
    time = []
    shear = []
    strain = []
    comp = []

    # READ FROM CSV FILE
    with open(filename, 'rb') as csvfile:
        for i in range(0,stVrst):
            next(csvfile)
        podatki = csv.reader(csvfile)
        for row in podatki:
            time.append(float(row[1]))
            shear.append(float(row[2]))
            strain.append(float(row[3]))
            comp.append(float(row[4]))

    # CONVERT LIST TO ARRAYS
    time = np.array(time)
    shear = np.array(shear)
    strain = np.array(strain)
    comp = np.array(comp)
    tc = time[0:len(time)/3]
    sc = strain[0:len(time)/3]
    cc =comp[0:len(time)/3]
    tr = time[len(time)/3:len(time)]
    sr = strain[len(time)/3:len(time)]

    return tc,tr,sc,sr,cc;

def burgerModel(params,t1,t2,y1=None,y2=None):
    l0   = params['l0'].value
    G0   = params['G0'].value
    l1   = params['l1'].value
    G1   = params['G1'].value

    tc   = 161
    tauC = params['tauC'].value

    strainC = tauC*t1/(l0*G0)+tauC/G0+tauC/G1*(1.0-np.exp(-t1/l1))

    gammaC = tauC*tc/(l0*G0)+tauC/G0+tauC/G1*(1.0-np.exp(-tc/l1))
    strainR = tauC*tc/(l0*G0)+tauC/G1*(np.exp(-(t2-tc)/l1))

    if (y1 is None) and (y2 is None):
        return np.r_[strainC,strainR]

    return np.r_[strainC-y1,strainR-y2];
# Burgersov model za voljnost
def burgerModelJ(params,t1,t2,y1=None,y2=None):
    l0   = params['l0'].value
    G0   = params['G0'].value
    l1   = params['l1'].value
    G1   = params['G1'].value

    J0 = 1.0/G0
    J1 = 1.0/G1
    tc   = 161
    tauC = params['tauC'].value

    JC = t1/(l0*G0)+J0+J1*(1.0-np.exp(-t1/l1))
    JR = tc/(l0*G0)+J1*(np.exp(-(t2-tc)/l1))

    if (y1 is None) and (y2 is None):
        return np.r_[JC,JR]

    return np.r_[JC-y1,JR-y2];

def KVModel(params,t1,t2,y1=None,y2=None):
    G = params['G'].value
    l = params['l'].value

    tauC = params['tauC'].value
    tc = 161

    strainC = tauC/G*(1.0-np.exp(-t1/l))

    gammaC = tauC/G*(1.0-np.exp(-tc/l))
    strainR = gammaC*(np.exp(-(t2-tc)/l))

    if (y1 is None) and (y2 is None):
        return strainC,strainR

    return np.r_[strainC-y1,strainR-y2];

def JeffreyModel(params,t1,t2,y1=None,y2=None):
    ni0   = params['ni0'].value
    l1   = params['l1'].value
    G1   = params['G1'].value

    tc   = 161
    tauC = params['tauC'].value

    strainC = tauC*t1/ni0+tauC/G1*(1.0-np.exp(-t1/l1))
    strainR = tauC*tc/ni0+tauC/G1*(np.exp(-(t2-tc)/l1))

    if (y1 is None) and (y2 is None):
        return np.r_[strainC,strainR]

    return np.r_[strainC-y1,strainR-y2];


def koefBurger(par):

    # EXTRACT PARAMETER VALUES
    p = par.valuesdict()
    koef = np.r_[p['G0'],p['l0'],p['G1'],p['l1']]

    return koef;

def koefKV(par):

    # EXTRACT PARAMETER VALUES
    p = par.valuesdict()
    koef = np.r_[p['G'],p['l']]

    return koef;

def koefJeffrey(par):

    # EXTRACT PARAMETER VALUES
    p = par.valuesdict()
    koef = np.r_[p['ni0'],p['G1'],p['l1']]

    return koef;
