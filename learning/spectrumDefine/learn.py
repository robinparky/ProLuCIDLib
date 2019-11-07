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

ACCURACY_THRESHOLD = 0.99

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
    with open (save + 'outputLabels', 'rb') as lp:
        outputLabels = np.array(pickle.load(lp))
    sp.close()
'''--------------------------------------------------------------------------'''

print(binArray.shape)
def generator(batchSize):
    while True:
        for i in range(0, len(binArray), batchSize):
            print(indexList[i:i + batchSize])
            yield binArray[i:i + batchSize], indexList[i:i + batchSize]



#implement callback function to stop training
# when accuracy reaches e.g. ACCURACY_THRESHOLD = 0.95
class myCallback(tf.keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs={}):
        if(logs.get('acc') > ACCURACY_THRESHOLD):
        	print("\nReached %2.2f%% accuracy, so stopping training!!" %(ACCURACY_THRESHOLD*100))
        	self.model.stop_training = True#Metrics for printing(mostly)


inputs = len(spectrums)
totalBins = len(binArray[0])
inputLayers = len(binArray)
outputLayers = len(outputLabels)


print("Inputs (spectrums): ",inputs)
print("\tShapes: ", binArray.shape);
print("\tType: ", binArray.dtype);
print("\tInput Layers: ", inputLayers, "\n")

print("Outputs:", )
print("\tNum: ", len(outputLabels));
print("\tType: ", outputLabels.dtype);
print("\tOutput Layers: ", outputLayers)

def create_model1():
    model = keras.Sequential()
    with tf.device('/gpu:0'):
        model.add(keras.layers.InputLayer(input_shape = (totalBins, )))
        model.add(keras.layers.Dense(outputLayers*.5, activation=tf.nn.relu))
    with tf.device('/gpu:1'):
        model.add(keras.layers.Dense(outputLayers*.5, activation=tf.nn.relu))
        model.add(keras.layers.Dense(outputLayers*.5, activation=tf.nn.relu))
    with tf.device('/gpu:2'):
        model.add(keras.layers.Dense(outputLayers*.5, activation=tf.nn.relu))
        model.add(keras.layers.Dense(outputLayers*.5, activation=tf.nn.relu))
    with tf.device('/gpu:3'):
        model.add(keras.layers.Dense(outputLayers, activation=tf.nn.softmax))
    return model



gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.333)

sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

def create_model2():
    model = keras.Sequential()
    model.add(keras.layers.InputLayer(input_shape = (totalBins, )))
    #model.add(keras.layers.Dense(outputLayers * 8, activation=tf.nn.relu))
    #model.add(keras.layers.Dense(outputLayers * 8, activation=tf.nn.relu))
    #model.add(keras.layers.Dense(outputLayers * .4, activation=tf.nn.relu))
    #model.add(keras.layers.Dense(outputLayers * .4, activation=tf.nn.relu))
    #model.add(keras.layers.Dense(outputLayers * .4, activation=tf.nn.relu))
    #model.add(keras.layers.Dense(outputLayers * .4, activation=tf.nn.relu))
    #model.add(keras.layers.Dense(outputLayers * .4, activation=tf.nn.relu))
    model.add(keras.layers.Dense(outputLayers * 2, activation=tf.nn.relu))
    model.add(keras.layers.Dense(outputLayers * 2, activation=tf.nn.relu))
    model.add(keras.layers.Dense(outputLayers * 2, activation=tf.nn.relu))
    model.add(keras.layers.Dense(outputLayers * 2, activation=tf.nn.relu))
    model.add(keras.layers.Dense(outputLayers, activation=tf.nn.softmax))

    return model



model  = create_model1()
model.summary()

#model = tf.keras.utils.multi_gpu_model(model, 4)
model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
'''--------------------------------------------------------------------------'''
#Train the data

#model.fit_generator(generator(batchSize), steps_per_epoch = np.ceil(len(binArray)/batchSize), epochs = numEpochs)
callbacks = myCallback()
model.fit(binArray, indexList, batch_size = batchSize, epochs = numEpochs, callbacks=[callbacks])
model.save_weights(outputPath)
#model.save(outputPath)


print ("Finished Training")
train = time.time()
print ("Elapsed Time: " + str(round(train - start)) + "\n")
