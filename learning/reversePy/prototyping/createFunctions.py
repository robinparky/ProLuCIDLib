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
    ions = {"b1": [], "b2": [], "bn1": [], "bn2": [], "bo1": [], "bo2": [], "y1": [], "y2": [], "yn1": []}
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

