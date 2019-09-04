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

attList = [];
original = []


labelList = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]);

for i, line in enumerate(file):
    originalPart = np.array(line.split("\t"))
    if i == 0:
        original = originalPart
    else:
        original = np.column_stack((original, originalPart))

    line = np.array(line.split("\t")[4:])
    if i == 0:
        attList = line
    else:
        attList = np.column_stack((attList, line))

#original = original.transpose()
#attList = attList.transpose()

print(original.shape)
print(attList.shape)


#data = pd.DataFrame(attList)
#print(data)

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

with open(outputPath +'testLabelList', 'wb') as sp:
    pickle.dump(testLabelList, sp, protocol=4)
sp.close()
writeData = time.time()
print ("Time for section: " + str(round(writeData - pullData)))
print ("Elapsed Time: " + str(round(writeData - start)) + "\n")
