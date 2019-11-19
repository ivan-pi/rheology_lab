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

    return time,shear,strain,comp;
