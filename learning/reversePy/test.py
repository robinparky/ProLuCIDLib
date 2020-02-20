import sqlite3
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd

import keras.backend as K
import sys
from keras.layers.convolutional import Conv1D
from keras.layers.core import Dense, Dropout, Masking
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model as keras_load_model
from keras.models import Sequential

def printResults(array1, array2):
    print("--------------------")
    for i in zip(array1, array2):
        print(i[0], " | ", i[1])


if len(sys.argv) != 2:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]


xTest = np.load(INPUT_PATH + "XTest.npy")
yTest = np.load(INPUT_PATH + "YTest.npy")
model = keras_load_model(INPUT_PATH + "model.h5")



predictions = model.predict(xTest)

print(yTest[0])

for i, val in enumerate(predictions):
    printResults(yTest[i], val)
