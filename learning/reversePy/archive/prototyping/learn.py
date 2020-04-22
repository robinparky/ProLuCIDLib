import sqlite3
import sys
import struct
import time
import pickle
import random
import math
import gc
import json
import pandas as pd
from pyteomics import mass as massC
import importlib
importlib.import_module("preprocessing")
importlib.import_module("options")
importlib.import_module("modeling")


from preprocessing import *
from options import *
from modeling import *

#from createFunctions.py import getIonMasses

import numpy as np


def split_train_validate(x, y, validate_percent=.33, seed=None):
    length = len(x)
    np.random.seed(seed)
    indexs = np.random.permutation(length)
    train_end = int((1 - validate_percent) * length)
    train_indexs = indexs[:train_end]
    validate_indexs = indexs[train_end:]
    x_train = x[train_indexs]
    y_train = y[train_indexs]
    x_validate = x[validate_indexs]
    y_validate = y[validate_indexs]
    return x_train, y_train, x_validate, y_validate, train_indexs, validate_indexs



def binarySearch (arr, l, r, x, tol):

    # Check base case
    if r >= l:
        mid = int(l + (r - l)/2)

        # If element is present at the middle itself
        if arr[mid][0] > x - tol and arr[mid][0] < x + tol:
            return mid

        # If element is smaller than mid, then it can only
        # be present in left subarray
        elif arr[mid][0] > x + tol:
            return binarySearch(arr, l, mid-1, x, tol)

        # Else the element can only be present in right subarray
        else:
            return binarySearch(arr, mid+1, r, x, tol)

    else:
        # Element is not present in the array
        return -1

#Healper Function, Converts bitstring to float
def convertFloat(element):
    return struct.iter_unpack('>f', element)

def lossConvert(loss, charge):
     if loss == '':
         return 0
     elif loss == 'n':
         return massC.calculate_mass(formula='NH3', charge = charge)
     elif loss == 'o':
         return massC.calculate_mass(formula='H2O', charge = charge)

def getIonMasses(peptide, types=('b', 'y'), maxcharge=2):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    ions = {"b1": [], "b2": [], "bn1": [], "bn2": [], "bo1": [], "bo2": [], "y1": [], "y2": [], "yn1": [], "yn2": [], "yo1": [], "yo2": [], }
    losses = ['', 'n', 'o']
    for ion_type in types:
        for charge in range(1, maxcharge+1):
            for lossT in losses:
                key = ion_type  + lossT + str(charge)
                loss = lossConvert(lossT, charge)
                for i in range(1, len(peptide)-1):
                    if ion_type[0] in 'abc':
                        ions[key].append( massC.fast_mass(peptide[:i], ion_type=ion_type, charge=charge))
                    else:
                        ions[key].append( massC.fast_mass(peptide[i:], ion_type=ion_type, charge=charge))
    return ions

def one_hot_encode(seq):
    codes = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    o = list(set(codes) - set(seq))
    s = pd.DataFrame(list(seq))
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    return a
    #show_matrix(a)
    e = a.values.flatten()
    return e

#Get time to calulate program time
start = time.time()

#Connect to database
conn = sqlite3.connect('db.db')

xtrainc2 = []
xtrainc3 = []
ytrainc2 = []
ytrainc3 = []

c2Arr = []
c3Arr = []

#peptideTemp = {"peptide":"", "modification": 'null', "ions":{}}
peptideTemp = {"peptide":"", "ions":{}}

peptideCnt = 0
specCnt = 0

print ("Pulling Data from Database")
c = conn.cursor()
c.execute("SELECT * FROM PeptideTable")
pepTable = c.fetchall()

cnt = 0

#Iterate through table and pull information about each peptide and its scans
#for ind in pepTable:
for j in range(15):

    #Match peptide in table to peptide in Spectra Table
    #pepID = (str(ind[0]), )
    pepID = (str(pepTable[j][0]), )

    #peptide = ind[1].split('.')[1]
    peptide = pepTable[j][1].split('.')[1]

    if '(' in peptide or 'Z' in peptide or 'B' in peptide or 'U' in peptide:
        continue

    #charge = int(ind[4])
    charge = pepTable[j][4]

    #print(pepID, ":" ,peptide)



    peptideCnt += 1;
    c.execute('SELECT *,rowid FROM SpectraTable WHERE peptideID=?', pepID)
    spectrums = c.fetchall()

    #Make sure peptide has more than 50 spectrums for accurate training.
    specCnt += 1 #increment peptide count

    #Search returned list of matched spectrums for each peptide
    for element in spectrums:
        ions = getIonMasses(peptide)

        spectrum = [] #Holder variable for spectrum

        #append the scan id
        #idList.append(element[4])
        #massList.append(float(element[3])/1000);
        #Grab mzArr and intArr from
        mzArr = list(convertFloat(element[1]))
        intArr = list(convertFloat(element[2]))

        for key, massList in ions.items():
            intensities = []
            for mass in massList:
                #print(key,":" ,mass)
                found = False
                mass = mass
                tolerance = 400 * mass/1000000
                #cdistance = 10000
                #cValue = 0

                result = binarySearch(mzArr, 0, len(mzArr) - 1, mass, tolerance)
                if result != -1:
                    intensities.append(intArr[result][0])
                else:
                    intensities.append(0)
            ions[key] = intensities


        tmp = peptideTemp
        tmp['peptide'] = peptide
        tmp['ions'] = ions
        #print(tmp)


        if charge == 2:
            #xtrainc2.append(peptide)
            #ytrainc2.append(ions)
            c2Arr.append(tmp)
        elif charge == 3:
            #xtrainc3.append(peptide)
            #ytrainc3.append(ions)
            c3Arr.append(tmp)
        else:
            cnt += 1
            #print("Charge:", charge, "invalid")



print(cnt, "Charges Invalid")

#xtrainc2 = np.array(xtrainc2)
#ytrainc2 = np.array(ytrainc2)
#xtrainc3 = np.array(xtrainc3)
#ytrainc3 = np.array(ytrainc3)

#c2 = pd.DataFrame({'peptide': xtrainc2, 'ions': ytrainc2})
#c3 = pd.DataFrame({'peptide': xtrainc3, 'ions': ytrainc3})

"""
print(c2.head())
print(c2[0]['peptide'])
for spec in c2:
    print(spec)
    #print(type(spec))
    print(spec['ions'])
"""

options = PeptideMS2Options()
cvt = DataConverter(options)

xc2, yc2 =cvt.data_to_tensor(c2Arr)
xc3, yc3 =cvt.data_to_tensor(c3Arr)

x_train, y_train, x_validate, y_validate, train_indexs, validate_indexs = split_train_validate(xc2, yc2)


model = build_model(options)
result = model.fit(x_train, y_train, epochs = 100, validation_data = (x_validate, y_validate))
print(result)

print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")
