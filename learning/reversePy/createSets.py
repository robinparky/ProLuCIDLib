import sqlite3
import sys
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd
import time

import keras.backend as K
from keras.layers.convolutional import Conv1D
from keras.layers.core import Dense, Dropout, Masking
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model as keras_load_model
from keras.models import Sequential



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
df = pd.read_pickle(INPUT_PATH)
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

max_len = 0
for i, pep in enumerate(modPeptides):
    if len(pep) > max_len:
        max_len = len(pep)

for i, pep in enumerate(modPeptides):
    while len(pep) < max_len:
        pep = pep + "0"
    modPeptides.iloc[i] = nlf_encode(pep)



df = pd.concat([df, modPeptides], axis = 1)
df.columns = ['sequence', 'ions', 'encodedSequence']
df = df[['sequence', 'encodedSequence', 'ions']]
df.head()
#print(np.array(df['encodedSequence'].iloc[0]).shape)


# In[9]:

print("Processing Ions")
ionsProcessed = df['ions'].copy()


maxCnt = 0
for i, ionDict in enumerate(ionsProcessed):
    for key,value in ionDict.items():
        if len(value) > maxCnt:
            maxCnt = len(value)
print(maxCnt)

for i, ionDict in enumerate(ionsProcessed):
    container  = []
    for key,value in ionDict.items():
        ionList = []
        for j in value:
            ionList.append(j)
        while len(ionList) < maxCnt:
            ionList.append(0)
        container.append(np.array(ionList))
    ionsProcessed.iloc[i] = np.array(container)
    #print(np.array(container).shape)
#ionsProcessed.head()

df = pd.concat([df, ionsProcessed], axis = 1)
df.columns = ['sequence', 'encodedSequence', 'ions', 'ionsProcessed']
print(df.iloc[0]['sequence'])
print(df.iloc[0]['encodedSequence'])
print(df.iloc[0]['ions'])
print(df.iloc[0]['ionsProcessed'])


print("Creating Training/Testing Sets")



inputArr = np.array(df['encodedSequence'])
inputArr = np.array([np.array(i) for i in inputArr])
outputArr = np.array(df['ionsProcessed'])
outputArr = np.array([np.array(i) for i in outputArr])
fragmentShape =  outputArr.shape[1]*outputArr.shape[2]

print("Sequence Array Shape: ", inputArr.shape)
print("ION Array  Shape: ", outputArr.shape)

inputArr = np.reshape(inputArr, (inputArr.shape[0], 1, inputArr.shape[1]))
outputArr = np.reshape(outputArr, (outputArr.shape[0],fragmentShape))

splitIdx = int(len(inputArr) * TRAIN_SET_PERCENTAGE)

xTrain = inputArr[:splitIdx]
xTest = inputArr[splitIdx:]

yTrain = outputArr[:splitIdx]
yTest = outputArr[splitIdx:]

print("XTrain Shape: ", xTrain.shape)
print("YTrain Shape: ", yTrain.shape)
print("XTest Shape: ", xTest.shape)
print("YTest Shape: ", yTest.shape)
print("\nDone ")
print("Time Taken: ", time.time() - start)

np.save(OUTPUT_PATH + "XTrain.npy", xTrain)
np.save(OUTPUT_PATH + "YTrain.npy", yTrain)
np.save(OUTPUT_PATH + "XTest.npy", xTest)
np.save(OUTPUT_PATH + "YTest.npy", yTest)
