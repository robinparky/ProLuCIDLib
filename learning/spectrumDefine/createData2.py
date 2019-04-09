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

"""
Bin constants
Define maximum y val before flooring
Define Number of y bins
Define the width of x bins

Calculate the width of a y bin.
"""


Y_TOP_BOUND  = 1 #Defined top bound
yBins = 6
BIN_WIDTH_Y = Y_TOP_BOUND/yBins
BIN_WIDTH_X = .5
#BIN_WIDTH_Y = 1

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
for ind in pepTable:
#for j in range(10,15):

    #Match peptide in table to peptide in Spectra Table
    peptide = (str(ind[0]), )
    #peptide = (str(pepTable[j][0]), )

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

print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")

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

#Calculate the number of xbins and ybins.
xBins = math.ceil(MAX_X / BIN_WIDTH_X)
yBins = math.ceil(Y_TOP_BOUND/ BIN_WIDTH_Y)
#yBins = math.ceil(MAX_Y / BIN_WIDTH_Y)

#Total number of bins
#totalBins = xBins + xBins * (yBins-1)
totalBins = xBins * yBins
print( "Putting all coords into ", totalBins, " bins.")
print( "XBINS: ", xBins, " YBINS: ", yBins)
print( "Bins are ", BIN_WIDTH_X, " by ", BIN_WIDTH_Y)
print( "Max X: ", MAX_X, " Max Y: ", MAX_Y, "\n")

#Initialize binarray and binList
binArray = []
binList = np.zeros(totalBins)

#Iterate through spectrumList and put every coordinate pair in a bin creating
# a new representation of the spectrum
for ind in range(0, len(spectrumList)):
    print(ind,"/",len(spectrumList), end="\r", flush = True)

    #Get number of points in the spectrum and iterate through each one
    numPoints = len(spectrumList[ind])
    for ele, i in enumerate(spectrumList[ind]):
        x = i[0] # X coord
        y = i[1] # Y coord

        #Floor to max val if necessary
        if  math.log10(y) >= Y_TOP_BOUND:
            y = Y_TOP_BOUND - 1
        else:
            y = math.floor(math.log10(y))
        #Calculate what bin to put coords in and increment one to the bin
        x = math.floor(x/BIN_WIDTH_X)
        #y = math.floor(y/BIN_WIDTH_Y)

        binNumber = int(x + (y*xBins))

        binList[binNumber] = 1 #CHANGED FROM +=, Uncomment for percentage bins
        if binList[binNumber] != 0:
            binList[binNumber] = (binList[binNumber] + (y/MAX_Y))/2
        else:
            binList[binNumber] = y/MAX_Y



    #Divide everything by the number of coords to get all bins to be [0,1]
    #binList = np.true_divide(binList,numPoints) #Uncomment for percentage bins

    """for i, ele in enumerate(binList):       #This is multiplying by intensity value
        if ele != 0:
            binList[i] = 10 ** (i//2000) * 5"""

    #Append the binList to a holder, reset the binlist to all 0's and loop
    binArray.append(binList)
    binList = np.zeros(totalBins)

#Reshape and print binarray
binArray = np.array(binArray)
binArray = np.reshape(binArray, (len(indexList), yBins, xBins, 1))
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

for i in range(0, len(testSpec)-1):
    if not np.array_equal(spectrumListSplit[0][i],testSpec[i]) or  not np.array_equal(binArraySplit[0][i],testBins[i]) or labelListSplit[0][i] != testLab[i] or indexListSplit[0][i] != testInd[i] or idListSplit[0][i] != testId[i]:
        print("---------------------------------")
        print("\nOriginal")
        print(spectrumListSplit[0][i][4])
        print(binArraySplit[0][i][4])
        print(labelListSplit[0][i])
        print(indexListSplit[0][i])
        print(idListSplit[0][i])
        print("\nTest")
        print(testSpec[i][4])
        print(testBins[i][4])
        print(testLab[i])
        print(testInd[i])
        print(testId[i])

for i in range(0, len(spectrumList)-1):
    if not np.array_equal(spectrumListSplit[1][i],spectrumList[i]) or not np.array_equal(binArraySplit[1][i],binArray[i]) or labelListSplit[1][i] != labelList[i] or indexListSplit[1][i] != indexList[i] or idListSplit[1][i] != idList[i]:
        print("---------------------------------")
        print("\nOriginal")
        print(spectrumListSplit[1][i][4])
        print(binArraySplit[1][i][4])
        print(labelListSplit[1][i])
        print(indexListSplit[1][i])
        print(idListSplit[1][i])
        print("\nTraining")
        print(spectrumList[i][4])
        print(binArray[i][4])
        print(labelList[i])
        print(indexList[i])
        print(idList[i])


#'''''''''''''''''''''''''''''''''''''''TEST Data
print("Peptides: ", peptideCnt, " spectrumList: ", len(spectrumList))
print("Saving Data")
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
