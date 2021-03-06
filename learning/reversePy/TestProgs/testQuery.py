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
                #for i in range(1, len(peptide)-1):
                for i in range(1, len(peptide)):
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
c.execute("select * from peptidetable")
peptable = c.fetchall()

cnt = 0
pepcnt = 0

starttime = time.time()

#iterate through table and pull information about each peptide and its scans
#for ind in peptable:
for j in range(1000):

    #match peptide in table to peptide in spectra table
    #pepid = (str(ind[0]), )
    pepid = (str(peptable[j][0]), )

    #peptide = ind[1].split('.')[1]
    peptide = peptable[j][1].split('.')[1]

    #peptidefull = ind[1]
    peptidefull = peptable[j][1]
    #print(len(peptide))

    if '(' in peptide or 'z' in peptide or 'b' in peptide or 'u' in peptide or 'x' in peptide:
        continue

    #charge = int(ind[4])
    charge = peptable[j][4]

    #print(pepid, ":" ,peptide)

    peptidecnt += 1;
    c.execute('select *,rowid from spectratable where peptideid=?', pepid)
    spectrums = c.fetchall()

    #make sure peptide has more than 50 spectrums for accurate training.
    speccnt += 1 #increment peptide count

    #search returned list of matched spectrums for each peptide
    for element in spectrums:
        ionlist = getionmasses(peptide)

        """
        if peptide == "sldgdlagr":
            print(ionlist)
            sys.exit()
        """

        ions = {}

        spectrum = [] #holder variable for spectrum

        #append the scan id
        #idlist.append(element[4])
        #masslist.append(float(element[3])/1000);
        #grab mzarr and intarr from
        mzarr = list(convertfloat(element[1]))
        intarr = list(convertfloat(element[2]))

        for key, masslist in ionlist.items():
            intensities = []
            for mass in masslist:
                #print(key,":" ,mass)
                found = false
                mass = mass
                tolerance = 400 * mass/1000000
                #cdistance = 10000
                #cvalue = 0

                result = binarysearch(mzarr, 0, len(mzarr) - 1, mass, tolerance)
                if result != -1:
                    intensities.append(intarr[result][0])
                else:
                    intensities.append(0)
            ions[key] = intensities

        tmp = peptidetemp
        tmp['peptidefull'] = peptidefull
        tmp['peptide'] = peptide
        tmp['masslist'] = ionlist
        tmp['ions'] = ions
        tmp['retentiontime'] = element[4]
        tmp['scan'] = element[6]
        tmp['filename'] = element[5]
        #print(tmp)


        if charge == 2:
            c2arr.append(copy.copy(tmp))
        elif charge == 3:
            c3arr.append(copy.copy(tmp))
        else:
            cnt += 1
            #print("charge:", charge, "invalid")
    pepcnt += 1
print("Done Binary Search Creation")
print("Time Taken: ", time.time() - startTime)


