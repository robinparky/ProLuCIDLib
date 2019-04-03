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
 with tf.device('/cpu:0'):
     with open (save + 'spectrumsList', 'rb') as sp:
         spectrums = pickle.load(sp)
     sp.close()
 '''--------------------------------------------------------------------------'''

 print("Training Data with ", numEpochs, " epochs.")
 print("Neural network is using ", batchSize, "as batchsize.
 print("\n")

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
