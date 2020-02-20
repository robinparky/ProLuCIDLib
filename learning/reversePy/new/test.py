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

def printResults(dictionary, array):
    print("--------------------")
    j = 0
    for key, value in dictionary.items():
        for i in value:
            print(array[j], " | ", i)
            j += 1


if len(sys.argv) != 2:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]


xTest = np.load(INPUT_PATH + "XTest.npy")
yTest = np.load(INPUT_PATH + "YTest.npy")
model = keras_load_model(INPUT_PATH + "model.h5")



predictions = model.predict(xTest)



for i, val in enumerate(predictions):
    dictionary = yTest[i]["ions"]
    printResults(dictionary, val)
