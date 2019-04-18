from __future__ import division
import sys
import pickle
import time
import random
import math
import gc


import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras

if len(sys.argv) != 4:
    print("Error with command line inputs")
    sys.exit(0)
else:
    inputPath = sys.argv[1]
    neuralPath = sys.argv[2]
    showResults = sys.argv[3]

'''--------------------------------------------------------------------------'''
print("Gathering Data from Files")
with open (inputPath + 'attList', 'rb') as sp:
    attList = pickle.load(sp)
sp.close()
with open (inputPath + 'labelList', 'rb') as sp:
    labelList = pickle.load(sp)
sp.close()
with open (inputPath + 'testAttList', 'rb') as sp:
    testAttList = pickle.load(sp)
sp.close()
with open (inputPath + 'testLabelList', 'rb') as sp:
    testLabelList = pickle.load(sp)
sp.close()
'''--------------------------------------------------------------------------'''

model = tf.keras.Sequential([
  tf.keras.layers.Dense(32, activation=tf.nn.relu, input_shape=(8,)),  # input shape required
  tf.keras.layers.BatchNormalization(),
  tf.keras.layers.Dense(16, activation=tf.nn.relu),
  tf.keras.layers.BatchNormalization(),
  tf.keras.layers.Dense(2, activation=tf.nn.sigmoid )
])

model.compile(optimizer='adam',
           loss='sparse_categorical_crossentropy',
           metrics=['accuracy'])

model.load_weights(neuralPath)
#model = tf.keras.models.load_model('saved/model.h5')

result = model.predict(testAttList)

#Test Acccuracy
np.array(testLabelList)
test_loss, test_acc = model.evaluate(testAttList, testLabelList)
print('Test accuracy:', test_acc)
