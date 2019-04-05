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
    inputPath = sys.argv[1] #Path to sqlite database
    outputPath = sys.argv[2] #Path to store created arrays
    testNumber = int(sys.argv[3]) #Number of test elements to create

#Get time to calulate program time
start = time.time()

file = open(inputPath, "r")

done = False
mBool = False # Boolean for in M bounds
normal = False


attList = []
labelList = []
data = []


maxIntensity = 0
minIntensity = 10000
maxNum = 0
minNum = 100000
maxLength = 0
minLength = 1000000



for line in file:
    line = line.split("\t")
    if line[0] == "S":
        chargeState = int(line[3])
        m1 = float(line[6])
        intensity = float(line[7])
        if intensity > maxIntensity:
            maxIntensity = intensity
        elif intensity < minIntensity:
            minIntensity = intensity
        pTrue = float(line[8])
        numMatched = int(line[9])
        if numMatched > maxNum:
            maxNum = numMatched
        elif numMatched < minNum:
            minNum = numMatched
        data.append(chargeState)
        #data.append(m1)
        data.append(intensity)
        data.append(pTrue)
        data.append(numMatched)
        done = False
    elif line[0] == "M" and mBool == False and done == False:
        m2 = float(line[3])
        xCorr = float(line[5])
        peptide = line[9]
        lenPeptide = len(peptide.split('.')[1])
        if lenPeptide > maxLength:
            maxLength = lenPeptide
        elif lenPeptide < minLength:
            minLength = lenPeptide
        normalized = abs((m1 - m2)/m1)
        normalized = math.modf(normalized)[0]

        #data.append(m2)
        data.append(xCorr)
        data.append(lenPeptide)
        data.append(normalized)
        mBool = True

    elif line[0] != "M" and mBool == True:
        if not (line[1][0] == "R" and line[1][1] == "e" and line[1][7] == "_"):
            normal = True
    elif line[0] == "M" and mBool == True:
        data.append(line[4])
        if normal == True:
            labelList.append((1, 0))
        else:
            labelList.append((0,1))
        #print(data)
        attList.append(data)

        data = []
        normal = False
        mBool = False
        done = True
print(maxIntensity)
print(minIntensity)
for i in range(len(attList)):
    attList[i][1] = (attList[i][1]-minIntensity)/(maxIntensity - minIntensity) * 1000
    attList[i][2]  = attList[i][2]
    attList[i][3] = (attList[i][3]- minNum)/(maxNum-minNum) * 100
    attList[i][4] = attList[i][4]
    attList[i][5] = (attList[i][5]-minLength)/(maxLength-minLength) * 10
    attList[i][6] = attList[i][6]
    attList[i][7] = float(attList[i][7]) * 10


for i in range(100):
    print(attList[i])
    print(labelList[i])

print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")



#Shuffle lists so that it is mised and same peptide specs are not adjacent
shuffle = list(zip(attList, labelList))
random.shuffle(shuffle)
attList, labelList = list(zip(*shuffle))
attList = list(attList)
labelList = list(labelList)

'''--------------------------------------------------------------------------'''

#Create test set
print("Creating Test Set")

attList = np.array(attList)
labelList = np.array(labelList)


#Split array into test and training sets
attListSplit = np.split(attList, [int(testNumber)])
labelListSplit = np.split(labelList, [int(testNumber)])

testAttList = attListSplit[0]
testLabelList = labelListSplit[0]

attList = attListSplit[1]
labelList = labelListSplit[1]


#'''''''''''''''''''''''''''''''''''''''Save Data


print("Created Test Database with ", len(attList), " elements")
print("Test size is of size ", len(testAttList))


with open(outputPath +'attList', 'wb') as sp:
    pickle.dump(attList, sp, protocol=4)
sp.close()

with open(outputPath +'labelList', 'wb') as sp:
    pickle.dump(labelList, sp, protocol=4)
sp.close()

with open(outputPath +'testAttList', 'wb') as sp:
    pickle.dump(testAttList, sp, protocol=4)
sp.close()

with open(outputPath +'testLabelList', 'wb') as sp:
    pickle.dump(testLabelList, sp, protocol=4)
sp.close()


writeData = time.time()
print ("Time for section: " + str(round(writeData - pullData)))
print ("Elapsed Time: " + str(round(writeData - start)) + "\n")
