import sys
import pickle

import sqlite3
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd

import keras.backend as K
from keras.layers.convolutional import Conv1D
from keras.layers.core import Dense, Dropout, Masking
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model as keras_load_model
from keras.models import Sequential

from keras.callbacks import Callback

if len(sys.argv) != 4 :
    print("Error with command line inputs")
    sys.exit(0)
else:
    databasePath = sys.argv[1] #Path to sqlite database
    outputPath = sys.argv[2] #Path to store created arrays
    testNumber = int(sys.argv[3]) #Number of test elements to create



class TerminateOnBaseline(Callback):
    """Callback that terminates training when either acc or val_acc reaches a specified baseline
    """
    def __init__(self, monitor='acc', baseline=0.9):
        super(TerminateOnBaseline, self).__init__()
        self.monitor = monitor
        self.baseline = baseline

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        acc = logs.get(self.monitor)
        if acc is not None:
            if acc >= self.baseline:
                print('Epoch %d: Reached baseline, terminating training' % (epoch))
                self.model.stop_training = True

# In[2]:


def binarySearch (arr, l, r, x, tol):

    # Check base case
    if r >= l:
        mid = int(l + (r - l)/2)

        # If element is present at the middle itself
        if arr[mid][0] > x - tol and arr[mid][0] < x + tol:
            return mid

        # If element is smaller than mid, then it can only
        # be present in left subarray
        elif arr[mid][0] > x + tol:
            return binarySearch(arr, l, mid-1, x, tol)

        # Else the element can only be present in right subarray
        else:
            return binarySearch(arr, mid+1, r, x, tol)

    else:
        # Element is not present in the array
        return -1


# In[3]:


#Healper Function, Converts bitstring to float
def convertFloat(element):
    return struct.iter_unpack('>f', element)


# In[4]:


def lossConvert(loss, charge):
     if loss == '':
         return 0
     elif loss == 'n':
         return massC.calculate_mass(formula='NH3', charge = charge)
     elif loss == 'o':
         return massC.calculate_mass(formula='H2O', charge = charge)

def getIonMasses(peptide, types=('b', 'y'), maxcharge=2):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    ions = {"b1": [], "b2": [], "bn1": [], "bn2": [], "bo1": [], "bo2": [], "y1": [], "y2": [], "yn1": [], "yn2": [], "yo1": [], "yo2": [], }
    losses = ['', 'n', 'o']
    for ion_type in types:
        for charge in range(1, maxcharge+1):
            for lossT in losses:
                key = ion_type  + lossT + str(charge)
                loss = lossConvert(lossT, charge)
                for i in range(1, len(peptide)-1):
                    if ion_type[0] in 'abc':
                        ions[key].append( massC.fast_mass(peptide[:i], ion_type=ion_type, charge=charge))
                    else:
                        ions[key].append( massC.fast_mass(peptide[i:], ion_type=ion_type, charge=charge))
    return ions


# In[5]:


#Connect to database
print(sys.argv)
conn = sqlite3.connect(databasePath)
xtrainc2 = []
xtrainc3 = []
ytrainc2 = []
ytrainc3 = []

c2Arr = []
c3Arr = []

#peptideTemp = {"peptide":"", "modification": 'null', "ions":{}}
peptideTemp = {"peptide":"", "ions":{}}

peptideCnt = 0
specCnt = 0

print ("Pulling Data from Database")
c = conn.cursor()
c.execute("SELECT * FROM PeptideTable")
pepTable = c.fetchall()

cnt = 0

#Iterate through table and pull information about each peptide and its scans
for ind in pepTable:
#for j in range(15):

    #Match peptide in table to peptide in Spectra Table
    pepID = (str(ind[0]), )
    #pepID = (str(pepTable[j][0]), )

    peptide = ind[1].split('.')[1]
    #peptide = pepTable[j][1].split('.')[1]

    if '(' in peptide or 'Z' in peptide or 'B' in peptide or 'U' in peptide:
        continue

    charge = int(ind[4])
    #charge = pepTable[j][4]

    #print(pepID, ":" ,peptide)

    peptideCnt += 1;
    c.execute('SELECT *,rowid FROM SpectraTable WHERE peptideID=?', pepID)
    spectrums = c.fetchall()

    #Make sure peptide has more than 50 spectrums for accurate training.
    specCnt += 1 #increment peptide count

    #Search returned list of matched spectrums for each peptide
    for element in spectrums:
        ions = getIonMasses(peptide)

        spectrum = [] #Holder variable for spectrum

        #append the scan id
        #idList.append(element[4])
        #massList.append(float(element[3])/1000);
        #Grab mzArr and intArr from
        mzArr = list(convertFloat(element[1]))
        intArr = list(convertFloat(element[2]))

        for key, massList in ions.items():
            intensities = []
            for mass in massList:
                #print(key,":" ,mass)
                found = False
                mass = mass
                tolerance = 400 * mass/1000000
                #cdistance = 10000
                #cValue = 0

                result = binarySearch(mzArr, 0, len(mzArr) - 1, mass, tolerance)
                if result != -1:
                    intensities.append(intArr[result][0])
                else:
                    intensities.append(0)
            ions[key] = intensities


        tmp = peptideTemp
        tmp['peptide'] = peptide
        tmp['ions'] = ions
        #print(tmp)


        if charge == 2:
            #xtrainc2.append(peptide)
            #ytrainc2.append(ions)
            c2Arr.append(tmp)
        elif charge == 3:
            #xtrainc3.append(peptide)
            #ytrainc3.append(ions)
            c3Arr.append(tmp)
        else:
            cnt += 1
            #print("Charge:", charge, "invalid")
print("Done")


# In[6]:


df = pd.DataFrame(c2Arr)
df = df.sample(frac=1)
df.head()


#
# #Encode Sequence
#
# #read the matrix a csv file on github
# nlf = pd.read_csv('https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/epitopepredict/mhcdata/NLF.csv',index_col=0)
#
# def show_matrix(m):
#     #display a matrix
#     cm = sns.light_palette("seagreen", as_cmap=True)
#     display(m.style.background_gradient(cmap=cm))
#
# def nlf_encode(seq):
#     x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)
#     #show_matrix(x)
#     e = x.values.flatten()
#     return e
#
# modPeptides = df['peptide'].copy()
# for i, pep in enumerate(modPeptides):
#     modPeptides.iloc[i] = nlf_encode(pep)
#
# """
# for i in modPeptides:
#     print(i.shape)
# """
#
# modPeptides.head()
# print(modPeptides.iloc[0])

# In[7]:

print("Encoding Sequences")
import seaborn as sns
codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#Encode Sequence
def show_matrix(m):
    #display a matrix
    cm = sns.light_palette("seagreen", as_cmap=True)
    display(m.style.background_gradient(cmap=cm))

def one_hot_encode(seq):
    o = list(set(codes) - set(seq))
    s = pd.DataFrame(list(seq))
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    #show_matrix(a)
    e = a.values.flatten()
    return e

modPeptides = df['peptide'].copy()
for i, pep in enumerate(modPeptides):
    modPeptides.iloc[i] = one_hot_encode(pep)

"""
for i in modPeptides:
    print(i.shape)
"""

modPeptides.head()
#print(modPeptides.iloc[0])


# In[8]:


df = pd.concat([df, modPeptides], axis = 1)
df.columns = ['sequence', 'ions', 'encodedSequence']
df = df[['sequence', 'encodedSequence', 'ions']]
df.head()
#print(np.array(df['encodedSequence'].iloc[0]).shape)


# In[9]:

print("Modifying Ions for input to Model")
ionsProcessed = df['ions'].copy()
for i, ionDict in enumerate(ionsProcessed):
    container  = []
    for key,value in ionDict.items():
        ionList = []
        for j in value:
            ionList.append(j)
        container.append(np.array(ionList))
    ionsProcessed.iloc[i] = np.array(container)
    #print(np.array(container).shape)
#ionsProcessed.head()

df = pd.concat([df, ionsProcessed], axis = 1)
df.columns = ['sequence', 'encodedSequence', 'ions', 'ionsProcessed']
df.head()


# In[13]:

print("Creating train/test sets")
inputArr = np.array(df['encodedSequence'])
inputArr = np.array([np.array(i) for i in inputArr])
outputArr = np.array([np.array(i) for i in df['ionsProcessed']])

inputArr = np.reshape(inputArr, (inputArr.shape[0], 1, inputArr.shape[1]))
outputArr = np.reshape(outputArr, (outputArr.shape[0], 240))

splitIndex = len(inputArr) - testNumber

xTrain = inputArr[:splitIndex]
xTest = inputArr[splitIndex:]

yTrain = outputArr[:splitIndex]
yTest = outputArr[splitIndex:]

print("xTrain: ", xTrain.shape)
print("yTrain: ", yTrain.shape)
print("xTest: ", xTest.shape)
print("yTest: ", yTest.shape)

with open(outputPath +'xTrain', 'wb') as sp:
    pickle.dump(xTrain, sp, protocol=4)
sp.close()

with open(outputPath +'xTest', 'wb') as sp:
    pickle.dump(xTest, sp, protocol=4)
sp.close()

with open(outputPath +'yTrain', 'wb') as sp:
    pickle.dump(yTrain, sp, protocol=4)
sp.close()

with open(outputPath +'yTest', 'wb') as sp:
    pickle.dump(yTest, sp, protocol=4)
sp.close()


# In[14]:

print("Creating and Printing Model")
model = Sequential()

model.add(Bidirectional(LSTM(3, input_shape = xTrain.shape)))
#model.add(Dense(20))
#model.add(Dropout(0.5))
#model.add(TimeDistributed(Dense(240, activation='relu')))
model.add(Dense(240, activation='relu'))

model.compile(
    loss="mean_squared_error",
    optimizer="adam")
model.fit(xTrain, yTrain, epochs = 10000, callbacks = [TerminateOnBaseline(monitor='acc', baseline=0.95)])
model.save(outputPath)
