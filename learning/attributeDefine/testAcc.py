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
with open (inputPath + 'spectrumsList', 'rb') as sp:
    spectrums = pickle.load(sp)
sp.close()
'''--------------------------------------------------------------------------'''

inputs = len(spectrums)
totalBins = len(binArray[0])
inputLayers = len(binArray)
output = np.unique(indexList)
outputLayers = len(set(indexList))

def create_model():
    model = keras.Sequential([
        keras.layers.InputLayer(input_shape = (totalBins, )),
        keras.layers.Dense(outputLayers * 8, activation=tf.nn.relu),
        keras.layers.Dense(outputLayers * 8, activation=tf.nn.relu),
        keras.layers.Dense(outputLayers, activation=tf.nn.softmax)
    ])
    return model
model = create_model()

model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])

model.load_weights(neuralPath)
#model = tf.keras.models.load_model('saved/model.h5')

result = model.predict(testBins)

#Test Acccuracy
np.array(testInd)
test_loss, test_acc = model.evaluate(testBins, testInd)
print('Test accuracy:', test_acc)
