import sqlite3
import numpy as np
import pandas as pd
import copy
import json
import sys
import time

import keras.backend as K
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

if len(sys.argv) != 4:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]
    MODEL_PATH = sys.argv[2]
    OUTPUT_PATH = sys.argv[3]

start = time.time()

peptideList = np.load(INPUT_PATH + "peptideList.npy")
predictArray = np.load(INPUT_PATH + "encodedPeptideList.npy")
#model = keras_load_model(INPUT_PATH + "model.h5", custom_objects={'cosine_similarity': cosine_similarity})
model = keras_load_model(MODEL_PATH)

predictions = model.predict(predictArray)

ionList = ["b1", "b2", "bn1","bn2", "bo1", "bo2", "y1", "y2", "yn1", "yn2", "yo1", "yo2" ]

jsonArray = []
tmpPred = {}

jsonArrayExp = []
tmpPredExp = {}


conn = sqlite3.connect(OUTPUT_PATH)
c = conn.cursor()
c.execute("CREATE TABLE IF NOT EXISTS predictions(Sequence TEXT,\
                                                  PrecursorMZ float,\
                                                  Retention_Time float,\
                                                  Charge integer,\
                                                  b1 TEXT,\
                                                  b2 TEXT,\
                                                  bn1 TEXT,\
                                                  bn2 TEXT,\
                                                  bo1 TEXT,\
                                                  bo2 TEXT,\
                                                  y1 TEXT,\
                                                  y2 TEXT,\
                                                  yn1 TEXT,\
                                                  yn2 TEXT,\
                                                  yo1 TEXT,\
                                                  yo2 TEXT\
                                                  )")

for i in range(len(predictArray)):
    peptide = peptideList[i]

    #Size of predictionArray split into the 12 sections
    split = len(predictions[0]) / 12

    #The size of the partition without padding.
    splitEnd = len(peptideList[0]) - 1

    ionObj = {}

    splitCnt = 0
    cnt = 0
    while splitCnt <12:
        done = False
        ionArray = []
        while True:
            if done == False:
                ionArray.append(str(predictions[i][cnt]))
            cnt += 1
            if cnt % split == 0:
                break
            elif cnt % split > splitEnd:
                done = True
        ionObj[ionList[splitCnt]] = ','.join(copy.copy(ionArray))
        splitCnt += 1

    c.execute("INSERT INTO predictions(Sequence, PrecursorMZ, Charge, Retention_Time,\
                                       b1,b2,bn1,bn2,bo1,bo2,y1,y2,yn1,yn2,yo1,yo2)\
                                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", \
                                       (peptide, 0, 2, 0, \
                                       ionObj["b1"], ionObj["b2"],\
                                       ionObj["bn1"], ionObj["bn2"],\
                                       ionObj["bo1"], ionObj["bo2"],\
                                       ionObj["y1"], ionObj["y2"],\
                                       ionObj["yn1"], ionObj["yn2"],\
                                       ionObj["yo1"], ionObj["yo2"]))
    conn.commit()

conn.close()
print("Done")
print("Elapsed Time: ", time.time() - start)
