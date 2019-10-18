import sqlite3
import sys
import struct
import time
import pickle
import random
import math
import gc
import matplotlib.pyplot as plt
import numpy as np

import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix

if len(sys.argv) != 4 :
    print("Error with command line inputs")
    sys.exit(0)
else:
    inputPath = sys.argv[1] #Path to sqlite database
    outputPath = sys.argv[2] #Path to store created arrays
    testNumber = int(sys.argv[3]) #Number of test elements to create

#Get time to calulate program time
start = time.time()

attList = []
labelList = []
data = []

base = "data/D2_pPro_B-0"
ext = ".sqt"

for i in range(1, 10):
    #file = open(inputPath, "r")
    file = open(base + str(i) + ext, "r")

    done = False
    mBool = False # Boolean for in M bounds
    normal = False

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

        elif line[0] != "M" and mBool == True and done == False:
            if not (line[1][0] == "R" and line[1][1] == "e" and line[1][7] == "_"):
                normal = True
        elif line[0] == "M" and mBool == True and done == False:
            data.append(line[4])
            if normal == True:
                labelList.append(1)
                #data.append(1)
            else:
                labelList.append(0)
                #data.append(0)
            attList.append(data)

            data = []
            normal = False
            mBool = False
            done = True

"""for i in range(len(attList)):
    attList[i][1] = (attList[i][1]-minIntensity)/(maxIntensity - minIntensity) * 1000
    attList[i][2]  = attList[i][2]
    attList[i][3] = (attList[i][3]- minNum)/(maxNum-minNum) * 100
    attList[i][4] = attList[i][4]
    attList[i][5] = (attList[i][5]-minLength)/(maxLength-minLength) * 10
    attList[i][6] = attList[i][6]
    attList[i][7] = float(attList[i][7]) * 10"""



print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")



#Shuffle lists so that it is mised and same peptide specs are not adjacent
shuffle = list(zip(attList, labelList))
random.shuffle(shuffle)
attList, labelList = list(zip(*shuffle))
attList = list(attList)
labelList = np.array(list(labelList))

attList = np.array(attList).astype(np.float)
attList.reshape((len(attList), 8))

print(attList.shape)
print(labelList.shape)

#attList = attList/attList.max(axis = 0)



#attList.sort(key=lambda x: x[4])
#for i, ele in enumerate(attList):
    #print(str(ele))



'''--------------------------------------------------------------------------'''
if 1 == 0:
    x1 = []
    x2 = []
    y11 = []
    y12 = []
    y21 = []
    y22 = []
    y31 = []
    y32 = []
    y41 = []
    y42 = []
    y51 = []
    y52 = []
    y61 = []
    y62 = []
    y71 = []
    y72 = []
    y81 = []
    y82 = []

    for i, ele in enumerate(labelList):
        if ele == 1:
            x1.append(1)
            y11.append(attList[i][0])
            y21.append(attList[i][1])
            y31.append(attList[i][2])
            y41.append(attList[i][3])
            y51.append(attList[i][4])
            y61.append(attList[i][5])
            y71.append(attList[i][6])
            y81.append(attList[i][7])
        if ele == 0:
            x2.append(0)
            y12.append(attList[i][0])
            y22.append(attList[i][1])
            y32.append(attList[i][2])
            y42.append(attList[i][3])
            y52.append(attList[i][4])
            y62.append(attList[i][5])
            y72.append(attList[i][6])
            y82.append(attList[i][7])

    plt.style.use('ggplot')
    plt.subplot(3,3,1)
    plt.title("Charge State")
    plt.scatter(x1,y11,c = 'r',s=2)
    plt.scatter(x2,y12,c = 'b', s=2)

    plt.subplot(3,3,2)
    plt.title("Intensity")
    plt.scatter(x1,y21,c = 'r',s=2)
    plt.scatter(x2,y22,c = 'b', s=2)

    plt.subplot(3,3,3)
    plt.title("[P True]")
    plt.scatter(x1,y31,c = 'r',s=2)
    plt.scatter(x2,y32,c = 'b', s=2)

    plt.subplot(3,3,4)
    plt.title("numMatched")
    plt.scatter(x1,y41,c = 'r',s=2)
    plt.scatter(x2,y42,c = 'b', s=2)

    plt.subplot(3,3,5)
    plt.title("xCorr")
    plt.scatter(x1,y51,c = 'r',s=2)
    plt.scatter(x2,y52,c = 'b', s=2)


    plt.subplot(3,3,6)
    plt.title("Length of Peptide")
    plt.scatter(x1,y61,c = 'r',s=2)
    plt.scatter(x2,y62,c = 'b', s=2)

    plt.subplot(3,3,7)
    plt.title("normalized")
    plt.scatter(x1,y71,c = 'r',s=2)
    plt.scatter(x2,y72,c = 'b', s=2)

    plt.subplot(3,3,8)
    plt.title("Delta CN")
    plt.scatter(x1,y81,c = 'r',s=2)
    plt.scatter(x2,y82,c = 'b', s=2)
    plt.show()

if 1 ==1:
    df = pd.DataFrame(attList, columns = ['Charge State', 'Intensity', 'P True', 'numMatched', 'xCorr', 'Length of Peptide', 'Normalized', 'DeltaCn'])
    scatter_matrix(data, alpha=0.2, figsize=(6, 6), diagonal='kde')

#Create test set
print("Creating Test Set")

attList = np.array(attList)
labelList = np.array(labelList)

"""
#Split array into test and training sets
attListSplit = np.split(attList, [int(testNumber)])
labelListSplit = np.split(labelList, [int(testNumber)])

testAttList = attListSplit[0]
testLabelList = labelListSplit[0]

attList = attListSplit[1]
labelList = labelListSplit[1]
"""
testAttList = attList
testLabelList = labelList


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
