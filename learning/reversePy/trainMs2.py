import sqlite3
from pyteomics import mass as massC
import struct
import sys
import numpy as np
import pandas as pd
import time

import keras.backend as K
from keras.layers.core import Dense, Dropout, Masking
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model as keras_load_model
from keras.models import Sequential
from keras.callbacks import EarlyStopping
from keras.utils import multi_gpu_model
start = time.time()

if len(sys.argv) != 4:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]
    BATCH_SIZE = int(sys.argv[2])
    EPOCHS = int(sys.argv[3])

def cosine_similarity(y_true, y_pred):
    length = K.int_shape(y_pred)[1]
    y_true = K.batch_flatten(y_true)
    y_pred = K.batch_flatten(y_pred)
    y_true = K.l2_normalize(y_true, axis=-1)
    y_pred = K.l2_normalize(y_pred, axis=-1)
    cos = K.sum(y_true * y_pred, axis=-1, keepdims=True)
    result = -K.repeat_elements(cos, rep=length, axis=1)
    return result

xTrain = np.load(INPUT_PATH + "XTrain.npy")
yTrain = np.load(INPUT_PATH + "YTrainMs2.npy")

fragmentShape =  yTrain.shape[1]

model = Sequential()
model.add(Bidirectional(LSTM(3, input_shape = xTrain.shape)))
model.add(Dense(300))
model.add(Dense(300))
model.add(Dropout(0.5))
model.add(Dense(300))
#model.add(TimeDistributed(Dense(fragmentShape, activation='relu')))
model.add(Dense(fragmentShape, activation='relu'))

#model = multi_gpu_model(model, gpus=4)

model.compile(
    loss="mean_squared_error",
    #optimizer="adam", metrics=[cosine_similarity])
    optimizer="adam")

model.fit(xTrain, yTrain,
          batch_size = BATCH_SIZE,
          epochs = EPOCHS)

print("Saving Model")
model.save(INPUT_PATH + "modelMs2.h5")

print("Done")
print("Total Time: ", time.time() - start)
