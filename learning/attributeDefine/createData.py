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
yBins = 1
BIN_WIDTH_Y = Y_TOP_BOUND/yBins
BIN_WIDTH_X = .5
#BIN_WIDTH_Y = 1


#Get time to calulate program time
start = time.time()

#Connect to database
conn = sqlite3.connect(databasePath)

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

#'''''''''''''''''''''''''''''''''''''''Save Data
print("Peptides: ", peptideCnt, " spectrumList: ", len(spectrumList))

with open(outputPath +'spectrumsList', 'wb') as sp:
    pickle.dump(spectrumList, sp, protocol=4)
sp.close()


writeData = time.time()
print ("Time for section: " + str(round(writeData - pullData)))
print ("Elapsed Time: " + str(round(writeData - start)) + "\n")
