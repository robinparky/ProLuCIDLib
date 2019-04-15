import sqlite3
import sys
import struct
import time
import pickle
import random
import math
import gc

import numpy as np

if len(sys.argv) != 4 :
    print("Error with command line inputs")
    sys.exit(0)
else:
    databasePath = sys.argv[1] #Path to sqlite database
    outputPath = sys.argv[2] #Path to store created arrays
    testNumber = int(sys.argv[3]) #Number of test elements to create


#Healper Function, Converts bitstring to float
def convertFloat(element):
    return struct.iter_unpack('>f', element)

#Helper function, sorts a list by x coordinate
def Sort(list):
    return(sorted(list, key = lambda x: x[1]))

#Get time to calulate program time
start = time.time()

#Connect to database
conn = sqlite3.connect(databasePath)

peptideCnt = 0 #Count of number of peptides for alter use
specCnt = 0
labelList = [] #List of label associated with peptide
spectrumList = [] #Spectrum to hold arrays of coords
idList = [] #List holds scan id associated with spectrum

print ("Pulling Data from Database")
c = conn.cursor()
c.execute("SELECT * FROM PeptideTable")
pepTable = c.fetchall()

#Iterate through table and pull information about each peptide and its scans
#for ind in pepTable:
for j in range(10,15):

    #Match peptide in table to peptide in Spectra Table
    #peptide = (str(ind[0]), )
    peptide = (str(pepTable[j][0]), )

    peptideCnt += 1;
    c.execute('SELECT *,rowid FROM SpectraTable WHERE peptideID=?', peptide)
    spectrums = c.fetchall()

    #Make sure peptide has more than 50 spectrums for accurate training.
    if len(spectrums) >= 50:
        specCnt += 1 #increment peptide count

        #Search returned list of matched spectrums for each peptide
        for element in spectrums:

            spectrum = [] #Holder variable for spectrum

            #append the scan id
            idList.append(element[4])

            #Grab mzArr and intArr from
            mzArr = convertFloat(element[1])
            intArr = convertFloat(element[2])

            #Iterate through both arrays and append coordinate pairs
            for i, j in zip(mzArr, intArr):
                coords = []
                mz = float(i[0])
                intensity = float(j[0])
                coords.append(mz)
                coords.append(intensity)
                spectrum.append(coords)

            #Sort the spectrum and append.
            spectrum = Sort(spectrum)
            spectrumList.append(spectrum)

            #Grab label corresponding peptide label from the peptide table and
            #append to the labelList
            c.execute('SELECT * FROM PeptideTable WHERE peptideID=?', peptide)
            peptideRow = c.fetchone()
            seq = peptideRow[14]
            labelList.append(seq)
    else:
        peptideCnt -= 1

#Convert list of labelList to list of indices. The indices will correspond to
# each unique peptide
indexList = []
noDuplicateLabels = list(set(labelList))
for i in noDuplicateLabels:
    for y in labelList:
        if i == y:
            indexList.append(noDuplicateLabels.index(i))
print(indexList)
print(len(set(labelList)))
print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")


'''--------------------------------------------------------------------------'''
#Shuffle lists so that it is mised and same peptide specs are not adjacent

shuffle = list(zip(spectrumList, labelList, indexList, idList))
random.shuffle(shuffle)
spectrumList, labelList, indexList, idList = list(zip(*shuffle))
spectrumList = list(spectrumList)
labelList = list(labelList)
indexList = list(indexList)
idList = list(idList)

'''--------------------------------------------------------------------------'''
#Find max x and y
MAX_X = 0
MAX_Y = 0
for j in spectrumList:
    for i in j:
        if i[0] > MAX_X:
            MAX_X = i[0]
        if i[1] > MAX_Y:
            MAX_Y = i[1]

#Bound each x and y to nearest 1000, ceiling to not cut off
MAX_X = int(math.ceil(MAX_X / 1000.0)) * 1000
#MAX_Y =int(math.ceil(MAX_Y / 1000.0)) * 1000
MAX_Y = math.ceil(MAX_Y)

binWidth = 2

xBins = int(MAX_X/binWidth)
print( "Max X: ", MAX_X, " Max Y: ", MAX_Y, "\n")
print("xBins: ",xBins)

#Initialize binarray and binList
binArray = []
binList = np.zeros(xBins)

#Iterate through spectrumList and put every coordinate pair in a bin creating
# a new representation of the spectrum
for ind in range(0, len(spectrumList)):
    print(ind,"/",len(spectrumList), end="\r", flush = True)

    #Get number of points in the spectrum and iterate through each one

    for i in spectrumList[ind]:
        x = i[0] # X coord
        y = i[1] # Y coord

        binNumber = int(math.floor(x/binWidth))
        if y > binList[binNumber]:
            binList[binNumber] = y

    #Append the binList to a holder, reset the binlist to all 0's and loop
    binList = np.true_divide(binList, MAX_Y)
    binArray.append(binList)
    binList = np.zeros(xBins)

#Reshape and print binarray
binArray = np.array(binArray)
binArray = np.reshape(binArray, (len(indexList), xBins))
print(binArray.shape)

#Create test set
print("Creating Test Set")

spectrumList = np.array(spectrumList)
labelList = np.array(labelList)
indexList = np.array(indexList)
idList = np.array(idList)


#Split array into test and training sets
spectrumListSplit = np.split(spectrumList, [int(testNumber)])
binArraySplit = np.split(binArray, [int(testNumber)])
labelListSplit = np.split(labelList, [int(testNumber)])
indexListSplit = np.split(indexList, [int(testNumber)])
idListSplit = np.split(idList, [int(testNumber)])

testSpec = spectrumListSplit[0]
testBins = binArraySplit[0]
testLab = labelListSplit[0]
testInd = indexListSplit[0]
testId = idListSplit[0]

spectrumList = spectrumListSplit[1]
binArray = binArraySplit[1]
labelList = labelListSplit[1]
indexList = indexListSplit[1]
idList = idListSplit[1]
#'''''''''''''''''''''''''''''''''''''''TEST Data
print("Peptides: ", peptideCnt, " spectrumList: ", len(spectrumList))

with open(outputPath +'spectrumsList', 'wb') as sp:
    pickle.dump(spectrumList, sp, protocol=4)
sp.close()

with open(outputPath +'labelList', 'wb') as sp:
    pickle.dump(labelList, sp, protocol=4)
sp.close()

with open(outputPath +'indexList', 'wb') as sp:
    pickle.dump(indexList, sp, protocol=4)
sp.close()

with open(outputPath +'idList', 'wb') as sp:
    pickle.dump(idList, sp, protocol=4)
sp.close()

with open(outputPath +'binArray', 'wb') as sp:
    pickle.dump(binArray, sp, protocol=4)
sp.close()
####################################################
with open(outputPath +'testSpec', 'wb') as sp:
    pickle.dump(testSpec, sp, protocol=4)
sp.close()

with open(outputPath +'testBins', 'wb') as sp:
    pickle.dump(testBins, sp, protocol=4)
sp.close()

with open(outputPath +'testInd', 'wb') as sp:
    pickle.dump(testInd, sp, protocol=4)
sp.close()

with open(outputPath +'testLab', 'wb') as sp:
    pickle.dump(testLab, sp, protocol=4)
sp.close()

with open(outputPath +'testId', 'wb') as sp:
    pickle.dump(testId, sp, protocol=4)
sp.close()

with open(outputPath +'outputLabels', 'wb') as sp:
    pickle.dump(noDuplicateLabels, sp, protocol=4)
sp.close()


writeData = time.time()
print ("Time for section: " + str(round(writeData - pullData)))
print ("Elapsed Time: " + str(round(writeData - start)) + "\n")
