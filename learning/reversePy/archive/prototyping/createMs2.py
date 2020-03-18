import sqlite3
import sys
import struct
import time
import pickle
import random
import math
import gc
import json
from pyteomics import mass as massC


import numpy as np

# Returns index of x in arr if present, else -1
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
                        ions[key].append( massC.fast_mass(peptide[:i], ion_type=ion_type, charge=charge) - loss)
                    else:
                        ions[key].append( massC.fast_mass(peptide[i:], ion_type=ion_type, charge=charge) - loss)
    return ions

#Get time to calulate program time
start = time.time()

#Connect to database
conn = sqlite3.connect('db.db')

container1 = []
container2 = []
container3 = []
container4 = []

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
for ind in pepTable:
#for j in range(15):

    #Match peptide in table to peptide in Spectra Table
    pepID = (str(ind[0]), )
    #pepID = (str(pepTable[j][0]), )

    peptide = ind[1].split('.')[1]
    #peptide = pepTable[j][1].split('.')[1]

    if '(' in peptide or 'Z' in peptide or 'B' in peptide or 'U' in peptide:
        continue

    charge = int(ind[4])
    #charge = pepTable[j][4]

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
            container2.append(tmp)
        elif charge == 3:
            container3.append(tmp)
        else:
            cnt += 1
            #print("Charge:", charge, "invalid")

print(cnt, "Charges Invalid")
"""
for i in peptideContainer:
    print("\n")
    print(i)
"""


with open('charge2.ions.json', 'w') as json_file:
      json.dump(container2, json_file)
with open('charge3.ions.json', 'w') as json_file:
      json.dump(container3, json_file)

print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")


