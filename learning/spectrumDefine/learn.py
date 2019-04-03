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

if len(sys.argv) != 5 :
    print(len(sys.argv))
    print("Error with command line inputs")
    sys.exit(0)
else:
    save = sys.argv[1]
    outputPath = sys.argv[2]
    batchSize = int(sys.argv[3])
    numEpochs = int(sys.argv[4])

start = time.time()

print("Running program")
print("Training Data with ", numEpochs, " epochs.")
print("\n")

'''--------------------------------------------------------------------------'''
print("Gathering Data from Files")
with tf.device('/cpu:0'):
    with open (save + 'spectrumsList', 'rb') as sp:
        spectrums = pickle.load(sp)
    sp.close()
    with open (save + 'labelList', 'rb') as lp:
        labelList = pickle.load(lp)
    sp.close()
    with open (save + 'indexList', 'rb') as lp:
        indexList = pickle.load(lp)
    sp.close()
    with open (save + 'idList', 'rb') as lp:
        idList = pickle.load(lp)
    sp.close()
    with open (save + 'binArray', 'rb') as lp:
        binArray = pickle.load(lp)
    sp.close()
'''--------------------------------------------------------------------------'''


#Metrics for printing(mostly)
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

"""    with tf.device('/cpu:0'):
        baseModel = keras.Sequential([
            keras.layers.InputLayer(input_shape = (totalBins, )),
            keras.layers.Dense(outputLayers * 4, activation=tf.nn.relu),
            keras.layers.Dense(outputLayers * 4, activation=tf.nn.relu),
            keras.layers.Dense(outputLayers, activation=tf.nn.softmax)
        ])
    model = keras.utils.multi_gpu_model(baseModel, gpus=2)"""



def create_model():
    baseModel = keras.Sequential()
    baseModel.add(keras.layers.InputLayer(input_shape = (totalBins, )))
    baseModel.add(keras.layers.Dense(outputLayers * 8, activation=tf.nn.relu))
    baseModel.add(keras.layers.Dense(outputLayers * 8, activation=tf.nn.relu))
    baseModel.add(keras.layers.Dense(outputLayers, activation=tf.nn.softmax))

    return baseModel

model  = create_model()
model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
#model.summary()
'''--------------------------------------------------------------------------'''



#Train the data
model.fit(binArray, indexList, batch_size = batchSize, epochs = numEpochs)

model.save_weights(outputPath)
#model.save(outputPath)


print ("Finished Training")
train = time.time()
print ("Elapsed Time: " + str(round(train - start)) + "\n")
