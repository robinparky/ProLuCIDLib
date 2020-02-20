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

nlf = pd.read_csv('https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/epitopepredict/mhcdata/NLF.csv',index_col=0)
def nlf_encode(seq):
    x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)
    #show_matrix(x)
    e = x.values.flatten()
    return e

encodedSequence = np.array(nlf_encode("GHEENGDAVTEPQVAEEK")).reshape(1,-1)

print(encodedSequence.shape)

predArr = []
predArr = np.array(encodedSequence)

print(predArr.shape)

predArr = np.reshape(predArr, (predArr.shape[0], 1, predArr.shape[1]))
"""
if len(sys.argv) != 2:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]


model = keras_load_model(INPUT_PATH + "model.h5")

predictions = model.predict(predArr)

print(predictions)
"""
