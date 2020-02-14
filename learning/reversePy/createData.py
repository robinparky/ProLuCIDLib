import sqlite3
import sys
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd
import time
import psutil
import copy


if len(sys.argv) != 3:
    print("Error with command line inputs")
    sys.exit(0)
else:
    DATABASE_PATH = sys.argv[1]
    OUTPUT_PATH = sys.argv[2]

#Binary Search for elements in mass Spectra within Certain Tolerance
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



#Connect to database
#conn = sqlite3.connect('/data/tyrande/data/MudPit/projects2012_01_05_20_21166.db')
#conn = sqlite3.connect('projects2012_01_05_20_21166.db')
#conn = sqlite3.connect('testLibDuplicateSpectra.db')
conn = sqlite3.connect(DATABASE_PATH)

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
pepCnt = 0

startTime = time.time()

#Iterate through table and pull information about each peptide and its scans
#for ind in pepTable:
for j in range(1000):

    #Match peptide in table to peptide in Spectra Table
    #pepID = (str(ind[0]), )
    pepID = (str(pepTable[j][0]), )

    #peptide = ind[1].split('.')[1]
    peptide = pepTable[j][1].split('.')[1]
    peptideFull = pepTable[j][1]
    #print(len(peptide))

    if '(' in peptide or 'Z' in peptide or 'B' in peptide or 'U' in peptide or 'X' in peptide:
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
        ionList = getIonMasses(peptide)
        ions = {}

        spectrum = [] #Holder variable for spectrum

        #append the scan id
        #idList.append(element[4])
        #massList.append(float(element[3])/1000);
        #Grab mzArr and intArr from
        mzArr = list(convertFloat(element[1]))
        intArr = list(convertFloat(element[2]))

        for key, massList in ionList.items():
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
        #tmp['peptideFull'] = peptideFull
        tmp['peptide'] = peptide
        #tmp['massList'] = ionList
        tmp['ions'] = ions
        #tmp['retentionTime'] = element[4]
        #tmp['scan'] = element[6]
        #tmp['fileName'] = element[5]
        #print(tmp)


        if charge == 2:
            c2Arr.append(copy.copy(tmp))
        elif charge == 3:
            c3Arr.append(copy.copy(tmp))
        else:
            cnt += 1
            #print("Charge:", charge, "invalid")
    pepCnt += 1
    if pepCnt % 1000 == 0:
        print(int(pepCnt /1000), "---", time.time() - startTime, "---", psutil.virtual_memory())
print("Done Pulling from Database, Saving Dataframes")

c2DF = pd.DataFrame(c2Arr)
c2DF = c2DF.sample(frac=1)
"""
for i in range(1000):
    obj = c2DF.iloc[i]
    print("----------------------------")
    print(i, obj.peptideFull)
    print("b1: ", obj.ions['b1'])
    print(obj['massList']['b1'])
    print("b2:" , obj.ions['b2'])
    print(obj['massList']['b2'])
    print("y1: ",obj.ions['y1'])
    print(obj['massList']['y1'])
    print("y2: ",obj.ions['y2'])
    print(obj['massList']['y2'])

    print("RT: ", obj.retentionTime)

    print("Scan: ", obj['scan'])
    print("FileName: ", obj['fileName'])
    print("\n")

"""
c3DF = pd.DataFrame(c3Arr)
c3DF = c2DF.sample(frac=1)

"""
for i in c2DF['peptide']:
    print(i)
"""
"""
for i in df["ions"]:
    for key, value in i.items():
        print(value)
"""

c2DF.to_pickle(OUTPUT_PATH + "c2DF.pkl")
c3DF.to_pickle(OUTPUT_PATH + "c3DF.pkl")
