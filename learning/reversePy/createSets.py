import sqlite3
import sys
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd
import time

from sklearn.preprocessing import minmax_scale

#Get command line inputs
if len(sys.argv) != 4:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]
    OUTPUT_PATH = sys.argv[2]
    TRAIN_SET_PERCENTAGE = float(sys.argv[3])


start = time.time()
print("Encoding Sequences")

#Encode Sequence
df = pd.read_pickle(INPUT_PATH).reindex()
#read the matrix a csv file on github
nlf = pd.read_csv('NLF.csv',index_col=0)

def show_matrix(m):
    #display a matrix
    cm = sns.light_palette("seagreen", as_cmap=True)
    display(m.style.background_gradient(cmap=cm))

def nlf_encode(seq):
    x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)
    #show_matrix(x)
    e = x.values.flatten()
    return e

modPeptides = df['peptide'].copy()
#print(modPeptides.head())

MAX_LEN = 50

keyDict = {}

for i, pep in enumerate(modPeptides):
    original = pep;
    #tmpCnt = max_len - len(pep)
    if pep not in keyDict:
        if len(pep) > MAX_LEN:
            df = df.drop[i]
            modPeptides = modPeptides.drop[i]
            continue
        while len(pep) < MAX_LEN:
            pep =pep + "0"
        encodedPeptide = nlf_encode(pep)
        modPeptides.iloc[i] = encodedPeptide
        keyDict[original] = encodedPeptide
    else:
        modPeptides.iloc[i] = keyDict[pep]


    """
    if i%100 == 0:
        print(i)
    """
#print(modPeptides.head())
#print(df.columns)

df = pd.concat([df, modPeptides], axis = 1)

#print(df.head())
df.columns = ['peptide', 'ions', 'peptideFull', 'massList', 'retentionTime', 'scan',
               'fileName', 'encodedSequence']
df = df[['peptide','encodedSequence', 'ions', 'peptideFull', 'massList', 'retentionTime', 'scan',
               'fileName']]
#print(df.head())
#print(np.array(df['encodedSequence'].iloc[0]).shape)


# In[9]:
print("Done")
print("Elapsed Time: ", time.time() - start)
print("Processing Ions")
ionsProcessed = df['ions'].copy()

"""
maxCnt = 0
for i, ionDict in enumerate(ionsProcessed):
    for key,value in ionDict.items():
        if len(value) > maxCnt:
            maxCnt = len(value)
"""
#print(maxCnt)

for i, ionDict in enumerate(ionsProcessed):
    container  = []
    for key,value in ionDict.items():
        ionList = []
        for j in value:
            ionList.append(j)
        while len(ionList) < MAX_LEN:
            ionList.append(0)
        #Normalize List here
        ionList = minmax_scale(ionList)
        container.append(np.array(ionList))
    ionsProcessed.iloc[i] = np.array(container)
    #print(np.array(container).shape)
#ionsProcessed.head()

df = pd.concat([df, ionsProcessed], axis = 1)

df.columns = ['peptide','encodedSequence', 'ions', 'peptideFull', 'massList', 'retentionTime', 'scan',
               'fileName', 'ionsProcessed']
df = df[['peptide','encodedSequence', 'ions','ionsProcessed', 'peptideFull', 'massList', 'retentionTime', 'scan',
               'fileName']]
"""
print("Printing")
print(df.iloc[len(df)-1]['peptide'])
print(df.iloc[len(df)-1]['encodedSequence'])
print(df.iloc[len(df)-1]['ions'])
print(df.iloc[len(df)-1]['ionsProcessed'])
"""

print("Done")
print("Elapsed Time: ", time.time() - start)
print("\n")

print("Creating Training/Testing Sets")

inputArr = np.array(df['encodedSequence'])
inputArr = np.array([np.array(i) for i in inputArr])

outputArrMs2 = np.array(df['ionsProcessed'])
outputArrMs2 = np.array([np.array(i) for i in outputArrMs2])

outputArrRt = np.array(df['retentionTime'])

fragmentShape =  outputArrMs2.shape[1]*outputArrMs2.shape[2]

print("Sequence Array Shape: ", inputArr.shape)
print("ION Array  Shape: ", outputArrMs2.shape)
print("Retention Array Shape: ", outputArrRt.shape)
print("\n")

inputArr = np.reshape(inputArr, (inputArr.shape[0], 1, inputArr.shape[1]))
outputArrMs2 = np.reshape(outputArrMs2, (outputArrMs2.shape[0],fragmentShape))

splitIdx = int(len(inputArr) * TRAIN_SET_PERCENTAGE)

xTrain = inputArr[:splitIdx]
xTest = inputArr[splitIdx:]

yTrainMs2 = outputArrMs2[:splitIdx]
yTestMs2 = outputArrMs2[splitIdx:]

yTrainRt = outputArrRt[:splitIdx]
yTestRt = outputArrRt[splitIdx:]


print("XTrain Shape: ", xTrain.shape)
print("YTrainMs2 Shape: ", yTrainMs2.shape)
print("YTrainRt Shape: ", yTrainRt.shape)
print("\n")

print("XTest Shape: ", xTest.shape)
print("YTestMs2 Shape: ", yTestMs2.shape)
print("YTestRt Shape: ", yTestRt.shape)

print("\nDone ")
print("Time Taken: ", time.time() - start)


dfTrain = df.iloc[:splitIdx]
dfTest = df.iloc[splitIdx:]

dfTrain.to_pickle(OUTPUT_PATH + "dfTrain.pkl")
dfTest.to_pickle(OUTPUT_PATH + "dfTest.pkl")

np.save(OUTPUT_PATH + "XTrain.npy", xTrain)
np.save(OUTPUT_PATH + "YTrainMs2.npy", yTrainMs2)
np.save(OUTPUT_PATH + "YTrainRt.npy", yTrainRt)

np.save(OUTPUT_PATH + "XTest.npy", xTest)
np.save(OUTPUT_PATH + "YTestMs2.npy", yTestMs2)
np.save(OUTPUT_PATH + "YTestRt.npy", yTestRt)
