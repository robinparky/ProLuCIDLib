from __future__ import division
import sys
import pickle
import time
import math

import gc



import tensorflow as tf
from tensorflow import keras
import numpy as np
import matplotlib.pyplot as plt

save = 'data/'

if len(sys.argv) == 2 :
    NUM_EPOCHS = int(sys.argv[1])
    neuralNetworkPath = 'saved/1dDenseWeights.h5'
elif(len(sys.argv) == 3):
    NUM_EPOCHS = int(sys.argv[1])
    neuralNetworkPath = sys.argv[2]
else:
    NUM_EPOCHS = 3
    neuralNetworkPath = 'saved/1dDenseWeights.h5'

start = time.time()

print("Running program")
print("Training Data with ", NUM_EPOCHS, " epochs.")
print("\n")

'''--------------------------------------------------------------------------'''
print("Gathering Data from Files")
with open (save + 'spectrumsList', 'rb') as sp:
    spectrums = pickle.load(sp)
sp.close()
with open (save + 'labelList', 'rb') as lp:
    labelList = pickle.load(lp)
sp.close()
with open (save + '/indexList', 'rb') as lp:
    indexList = pickle.load(lp)
sp.close()
with open (save + '/idList', 'rb') as lp:
    idList = pickle.load(lp)
sp.close()
with open (save + '/binArray', 'rb') as lp:
    binArray = pickle.load(lp)
sp.close()
'''--------------------------------------------------------------------------'''

inputs = len(spectrums)
totalBins = len(binArray[0])
inputLayers = len(binArray)
output = np.unique(indexList)
outputLayers = len(set(indexList))


print("Inputs (spectrums): ",inputs)
print("\tShapes: ", binArray.shape);
print("\tType: ", binArray.dtype);
print("\tInput Layers: ", inputLayers, "\n")

print("Outputs:", )
print("\tShapes: ", output.shape);
print("\tType: ", output.dtype);
print("\tOutput Layers: ", outputLayers)

"""print("Training Inputs: ")
print("\tShapes: ", testBins.shape);
print("\tType: ", testBins.dtype);
print("\tInput Layers: ", inputLayers, "\n")
print("Outputs:", )
print("\tShapes: ", testInd.shape);
print("\tType: ", testInd.dtype);
print("\tOutput Layers: ", outputLayers)
print("\n")"""

"""def create_model():

    with tf.device('/cpu:0'):
        baseModel = keras.Sequential([
            keras.layers.InputLayer(input_shape = (totalBins, )),
            keras.layers.Dense(inputLayers/2, activation=tf.nn.relu),
            #keras.layers.Dense(inputLayers/2, activation=tf.nn.relu),
            keras.layers.Dense(outputLayers, activation=tf.nn.softmax)
        ])
    model = keras.utils.multi_gpu_model(baseModel, gpus=2)
    return model"""

def create_model():

    with tf.device('/cpu:0'):
        baseModel = keras.Sequential([
            keras.layers.InputLayer(input_shape = (totalBins, )),
            keras.layers.Dense(inputLayers/2, activation=tf.nn.relu),
            keras.layers.Dense(inputLayers/2, activation=tf.nn.relu),
            keras.layers.Dense(outputLayers, activation=tf.nn.softmax)
        ])
    model = keras.utils.multi_gpu_model(baseModel, gpus=2)
    return model

model  = create_model()
model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
model.summary()
'''--------------------------------------------------------------------------'''

#Train the data
model.fit(binArray, indexList, batch_size = 32, epochs = NUM_EPOCHS)

model.save_weights(neuralNetworkPath)


print ("Finished Training")
train = time.time()
print ("Time for section: " + str(round(train - start)))
print ("Elapsed Time: " + str(round(train - start)) + "\n")