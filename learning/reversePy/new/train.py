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



if len(sys.argv) != 4:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]
    BATCH_SIZE = sys.argv[2]
    EPOCHS = sys.argv[3]

xTrain = np.load(INPUT_PATH + "XTrain.npy")
yTrain = np.load(INPUT_PATH + "YTrain.npy")

model = Sequential()

model.add(Bidirectional(LSTM(3, input_shape = xTrain.shape)))
model.add(Dense(100))
model.add(Dropout(0.5))
#model.add(TimeDistributed(Dense(fragmentShape, activation='relu')))
model.add(Dense(fragmentShape, activation='relu'))

model.compile(
    loss="mean_squared_error",
    optimizer="adam")
model.fit(xTrain, yTrain, batch_size = BATCH_SIZE, epochs = EPOCHS)

model.save(INPUT_PATH + "model.h5")
