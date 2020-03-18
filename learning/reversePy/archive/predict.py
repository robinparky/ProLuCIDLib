import sqlite3
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd
import copy
import json

import keras.backend as K
import sys
from keras.layers.convolutional import Conv1D
from keras.layers.core import Dense, Dropout, Masking
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model as keras_load_model
from keras.models import Sequential

def cosine_similarity(y_true, y_pred):
    length = K.int_shape(y_pred)[1]
    y_true = K.batch_flatten(y_true)
    y_pred = K.batch_flatten(y_pred)
    y_true = K.l2_normalize(y_true, axis=-1)
    y_pred = K.l2_normalize(y_pred, axis=-1)
    cos = K.sum(y_true * y_pred, axis=-1, keepdims=True)
    result = -K.repeat_elements(cos, rep=length, axis=1)
    return result

def printResults(array1, array2):
    print("--------------------")
    for i in zip(array1, array2):
        print(i[0], " | ", i[1])

if len(sys.argv) != 2:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]

df = pd.read_pickle(INPUT_PATH + "dfTest.pkl").reindex()

xTest = np.load(INPUT_PATH + "XTest.npy")
yTest = np.load(INPUT_PATH + "YTest.npy")
#model = keras_load_model(INPUT_PATH + "model.h5", custom_objects={'cosine_similarity': cosine_similarity})
model = keras_load_model(INPUT_PATH + "model.h5")



predictions = model.predict(xTest)

ionList = ["b1", "b2", "bn1","bn2", "bo1", "bo2", "y1", "y2", "yn1", "yn2", "yo1", "yo2" ]

"""
for i, val in enumerate(predictions):
    printResults(yTest[i], val)
"""

print(xTest.shape)
print(yTest.shape)
print(df.shape)

jsonArray = []
tmpPred = {}

jsonArrayExp = []
tmpPredExp = {}

for i in range(len(yTest)):
    target = df.iloc[i]
    print(target["peptideFull"], target["scan"], target["fileName"])
    print(target["ions"])
    split = len(yTest[0]) / 12
    splitEnd = len(target["ions"]["b1"]) - 1
    #print(split)

    tmpPred['peptide'] = target["peptideFull"].split('.')[1]

    ionObj = {}
    ionObjExp = {}

    splitCnt = 0
    cnt = 0
    while splitCnt <12:
        print(ionList[splitCnt], "----------------")
        done = False
        ionArray = []
        ionArrayExp = []
        while True:
            if done == False:
                print(yTest[i][cnt], "\t" ,predictions[i][cnt])

                ionArray.append(float(predictions[i][cnt]))
                ionArrayExp.append(float(yTest[i][cnt]))
            cnt += 1
            if cnt % split == 0:
                break
            elif cnt % split > splitEnd:
                done = True
        ionObj[ionList[splitCnt]] = copy.copy(ionArray)
        ionObjExp[ionList[splitCnt]] = copy.copy(ionArrayExp)

        splitCnt += 1
    tmpPred['ions'] = copy.copy(ionObj)
    tmpPredExp['ions'] = copy.copy(ionObjExp)
    jsonArray.append(copy.copy(tmpPred))
    jsonArrayExp.append(copy.copy(tmpPredExp))

print("Done")

print(jsonArray[0], jsonArray[1])
print(jsonArrayExp[0], jsonArrayExp[1])

with open('predictions.json', 'w', encoding='utf-8') as f:
    json.dump(jsonArray, f, ensure_ascii=False, indent=4)

with open('experimental.json', 'w', encoding='utf-8') as g:
    json.dump(jsonArrayExp, g, ensure_ascii=False, indent=4)
