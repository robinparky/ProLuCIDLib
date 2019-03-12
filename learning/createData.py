import sqlite3
import sys
import struct
import time
import pickle
import random
import math
import gc

import numpy as np

#CONSTANTS
NUMBER_PEPTIDES = 30
NUMBER_TO_TEST = 100

Y_TOP_BOUND  = 100000
yBins = 10

BIN_WIDTH_Y = Y_TOP_BOUND/yBins
#BIN_WIDTH_Y = 10
BIN_WIDTH_X = 1

def convertFloat(element):
    return struct.iter_unpack('>f', element)
def Sort(list):# Sorts by second element in sublist
    return(sorted(list, key = lambda x: x[1]))
if len(sys.argv) == 2 :
    databasePath = sys.argv[1]
    outputPath = "data/"
elif(len(sys.argv) == 3):
    databasePath = sys.argv[1]
    outputPath = sys.argv[2]
else:
    databasePath = "data/testLibDuplicateSpectra.db"
    outputPath = "data/"

start = time.time()
conn = sqlite3.connect(databasePath)

zeroFill = 0;
specZeroFill = 0;
cnt = 0
c = conn.cursor()
print ("Pulling Data from Database")
c.execute("SELECT * FROM PeptideTable")
pepTable = c.fetchall()
labelList = []
numLabels = []
spectrumList = []
idList = []
#for ind in pepTable:
for j in range(0, NUMBER_PEPTIDES):
    #peptide = (str(ind[0]), )
    peptide = (str(pepTable[j][0]), ) #Match peptide in table to peptide in Spectra Table
    c.execute('SELECT *,rowid FROM SpectraTable WHERE peptideID=?', peptide)
    spectrums = c.fetchall()
    if len(spectrums) >= 50:
        cnt += 1
        #print(len(spectrums))
        for element in spectrums: #Search returned list of matched spectrums
            idList.append(element[4])
            spectrum = [];
            mzArr = convertFloat(element[1])
            intArr = convertFloat(element[2])
            for i, j in zip(mzArr, intArr):
                coords = []
                mz = float(i[0])
                intensity = float(j[0]);
                #print(mz, " | ", intensity)
                coords.append(mz);
                coords.append(intensity);
                spectrum.append(coords);

            spectrum = Sort(spectrum)
            spectrumList.append(spectrum)
            """for i in range(0, len(spectrums)):
                spectrumList.append(spectrum);"""

            c.execute('SELECT * FROM PeptideTable WHERE peptideID=?', peptide)
            peptideRow = c.fetchone()
            seq = peptideRow[14]
            labelList.append(seq);
            """for i in range(0, len(spectrums)):
                labelList.append(seq);
            break"""
            #Print Data/

            """print("Protein: ", seq);
            print("Number Coords: " + str(len(spectrum)))
            print("Spectrum: ")
            for i in spectrum:
                line = '{:<24}  {:<24}'.format(str(i[0]), str(i[1]))
                print(line)"""

#Convert list of labelList to list of numberss
indexList = []
noDuplicateLabels = list(set(labelList));
for i in noDuplicateLabels:
    for y in labelList:
        if i == y:
            indexList.append(noDuplicateLabels.index(i))



print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")

#Shuffle Lists
shuffle = list(zip(spectrumList, labelList, indexList, idList))
random.shuffle(shuffle);
spectrumList, labelList, indexList, idList = list(zip(*shuffle));
spectrumList = list(spectrumList);
labelList = list(labelList);
indexList = list(indexList);
idList = list(idList);

'''--------------------------------------------------------------------------'''
#Find max x and y1
MAX_X = 0;
MAX_Y = 0;
for j in spectrumList:
    for i in j:
        if i[0] > MAX_X:
            MAX_X = i[0]
        if i[1] > MAX_Y:
            MAX_Y = i[1]


MAX_X = int(math.ceil(MAX_X / 1000.0)) * 1000
MAX_Y =int(math.ceil(MAX_Y / 1000.0)) * 1000

xBins = math.ceil(MAX_X / BIN_WIDTH_X);
yBins = math.ceil(Y_TOP_BOUND/ BIN_WIDTH_Y);
#yBins = math.ceil(MAX_Y / BIN_WIDTH_Y);

totalBins = xBins + xBins * (yBins-1)
print( "Putting all coords into ", totalBins, " bins.")
print( "XBINS: ", xBins, " YBINS: ", yBins)
print( "Bins are ", BIN_WIDTH_X, " by ", BIN_WIDTH_Y)
print( "Max X: ", MAX_X, " Max Y: ", MAX_Y, "\n")
binArray = []
binList = np.zeros(totalBins);
yBinTemp = []

for ind in range(0, len(spectrumList)):
    print(ind,"/",len(spectrumList), end="\r", flush = True)
    numPoints = len(spectrumList[ind])
    for cnt, i in enumerate(spectrumList[ind]):
        """print("-----------------------------------------")
        print("MAX_X: ", MAX_X, " MAX_Y: ", MAX_Y)
        print("Number Points: ", numPoints)
        print("BINS: ", xBins, " ", yBins)
        print ("BIN Width: ", BIN_WIDTH_X, BIN_WIDTH_Y)
        print("Binlength: ", len(binList))
        print("x val: ", i[0])
        print("y val: ", i[1])"""
        if i[0] == 0 and i[1] == 0:
            numPoints -= 1
        else:
            x = i[0]
            y = i[1]
            if  y > Y_TOP_BOUND:
                y = Y_TOP_BOUND - 1
            x = math.floor(x/BIN_WIDTH_X);
            y = math.floor(y/BIN_WIDTH_Y);from pprint import pprint
            binNumber = int(x + (y*xBins));
            #print(x, " | ", y)
            #print(binNumber)
            binList[binNumber] += 1
    binList = np.true_divide(binList,numPoints)
    binArray.append(binList)
    binList = np.zeros(totalBins)
binArray = np.array(binArray)
binArray = np.reshape(binArray, (len(indexList), totalBins))
print(binArray.shape)
testSpec = []
testBins = []
testInd = []
testLab = []
testId = []
gc.collect()
#Create test set
spectrumList = np.array(spectrumList)
labelList = np.array(labelList)
indexList = np.array(indexList)
idList = np.array(idList)
print("Creating Test Set")
for i in range(0, NUMBER_TO_TEST):
    print(i,"/",NUMBER_TO_TEST, end="\r", flush = True)
    testSpec.append(spectrumList[0])
    testBins.append(binArray[0])
    testInd.append(indexList[0])
    testLab.append(labelList[0])
    testId.append(idList[0])
    np.delete(spectrumList, 0)
    np.delete(binArray, 0)
    np.delete(indexList, 0)
    np.delete(labelList, 0)
    np.delete(idList, 0)

testSpec = np.array(testSpec)
testBins = np.array(testBins)
testInd = np.array(testInd)
testLab = np.array(testLab)
testId = np.array(testId)

#'''''''''''''''''''''''''''''''''''''''TEST Data
print("Peptides: ", cnt, " spectrumList: ", len(spectrumList))
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


writeData = time.time()
print ("Time for section: " + str(round(writeData - pullData)))
print ("Elapsed Time: " + str(round(writeData - start)) + "\n")
