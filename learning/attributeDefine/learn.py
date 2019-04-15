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
print("\n")

'''--------------------------------------------------------------------------'''
print("Gathering Data from Files")
with open (save + 'attList', 'rb') as sp:
    attList = pickle.load(sp)
sp.close()
with open (save + 'labelList', 'rb') as sp:
    labelList = pickle.load(sp)
sp.close()
with open (save + 'testAttList', 'rb') as sp:
    testAttList = pickle.load(sp)
sp.close()
with open (save + 'testLabelList', 'rb') as sp:
    testLabelList = pickle.load(sp)
sp.close()
'''--------------------------------------------------------------------------'''

print("Training Data with ", numEpochs, " epochs.")
print("Neural network is using ", batchSize, "as batchsize.")
print("\n")

inputs = len(attList[0])

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
model.summary()
'''--------------------------------------------------------------------------'''
print(np.shape(attList))
print(np.shape(labelList))

#Train the data
model.fit(attList, labelList, batch_size = batchSize, epochs = numEpochs)

model.save_weights(outputPath)
#model.save(outputPath)


print ("Finished Training")
train = time.time()
print ("Elapsed Time: " + str(round(train - start)) + "\n")
