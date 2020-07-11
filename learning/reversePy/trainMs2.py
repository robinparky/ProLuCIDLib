import sqlite3
from pyteomics import mass as massC
import struct
import sys
import numpy as np
import pandas as pd
import time

import keras.backend as K
from keras.layers.core import Dense, Masking
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed

from keras.models import load_model as keras_load_model
from keras.models import Sequential
from keras.utils import multi_gpu_model
start = time.time()

if len(sys.argv) != 4:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]
    BATCH_SIZE = int(sys.argv[2])
    EPOCHS = int(sys.argv[3])

xTrain = np.load(INPUT_PATH + "XTrain.npy")
yTrain = np.load(INPUT_PATH + "YTrainMs2.npy")

fragmentShape =  yTrain.shape[1]

model = Sequential()
model.add(Bidirectional(LSTM(3, input_shape = xTrain.shape)))
model.add(Dense(300))
model.add(Dense(300))
model.add(Dropout(0.5))
model.add(Dense(300))
model.add(Dense(fragmentShape, activation='relu'))

#model = multi_gpu_model(model, gpus=4)

model.compile(
    loss="mean_squared_error",
    optimizer="adam")

model.fit(xTrain, yTrain,
          batch_size = BATCH_SIZE,
          epochs = EPOCHS)

print("Saving Model")
model.save(INPUT_PATH + "modelMs2.h5")

print("Done")
print("Total Time: ", time.time() - start)
