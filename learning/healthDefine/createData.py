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
    inputPath = sys.argv[1] #Path to text data file
    outputPath = sys.argv[2] #Path to store created arrays
    testNumber = int(sys.argv[3]) #Number of test elements to create

#Get time to calulate program time
start = time.time()

file = open(inputPath, "r")

labelList = [];
attList = [];

for i, line in enumerate(file):
    if i != 0:
        tmp = []
        line = line.split("\t")
        labelList.append(int(line[0]))
        for ele in range(1, len(line)):
            tmp.append(float(line[ele]))
        attList.append(tmp)
"""
print(labelList)
for i, ele in enumerate(attList):
    print(i)
    print(str(ele[0]), labelList[i])
"""

attList = np.array(attList)
print(attList.shape)

print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")

si = []

"""
for i in range(len(attList[0])):
    sum = 0
    #print(i)
    for g in range(36):
        xprime = 0
        if i == 0:
            xprime = attList[g][i+1]
        elif i == len(attList[0]) - 1:
            xprime =  attList[g][i-1]
        else:
            if(abs(attList[g][i] - attList[g][i+1]) > abs(attList[g][i] - attList[g][i-1])):
                xprime = attList[g][i+1]
            else:
                xprime = attList[g][i-1]
        sum += (attList[g][i] + xprime + 1)%2
    sum = sum/36
    si.append(sum)
for i in si:
    print(i)
"""

for i in range(len(attList[0])):
    x = []
    y = []
    for g in range(36):
        print(attList[g][i])
        if labelList[g] == 1:
            y.append(0)
            x.append(attList[g][i])
        elif labelList[g] == 0:
            y.append(1)
            x.append(attList[g][i])

    #x.sort()
    sorted = list(zip(x,y))
    random.shuffle(sorted)
    x,y = list(zip(*sorted))
    x = list(x)
    y = list(y)
    print(x)
    plt.xticks(rotation='vertical')
    plt.scatter(x,y)
    plt.show()



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


"""
testAttList = attList
testLabelList = labelList
"""

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
