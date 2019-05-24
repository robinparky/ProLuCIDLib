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
    tf.keras.layers.Dense(len(attList[0]), activation=tf.nn.relu, input_shape=(len(attList[0]),)),  # input shape required
    tf.keras.layers.Dense(300, activation=tf.nn.relu),
    tf.keras.layers.Dense(200, activation=tf.nn.relu),
    tf.keras.layers.Dense(100, activation=tf.nn.relu),
    tf.keras.layers.Dense(50, activation=tf.nn.relu),
    tf.keras.layers.Dense(25, activation=tf.nn.relu),
    tf.keras.layers.Dense(1, activation=tf.nn.sigmoid )
])
model.compile(optimizer='adam',
           loss='sparse_categorical_crossentropy',
           metrics=['accuracy'])

model.load_weights(neuralPath)
#model = tf.keras.models.load_model('saved/model.h5')

result = model.predict(testAttList)

outputLabels =  [0,1]

correctCnt = 0
cnt = 0

resultList = []

print('\n')
print("1 Denotes Sick, 0 Denotes Healthy")
print('\n')
print("  Predicted  |  Actual ")

for i, ele in enumerate(result):
    answer = 0
    if ele[0] > .5:
        answer = 1
    print(i, "      ",answer, "  |   ", testLabelList[i])
    if answer == testLabelList[i]:
        correctCnt += 1
    cnt += 1

print('\n')
print("Tested", correctCnt, "correct out of", cnt)
print("Accuracy:", str(correctCnt/cnt))

#Test Acccuracy
"""
np.array(testLabelList)
test_loss, test_acc = model.evaluate(testAttList, testLabelList)
print('Test accuracy:', test_acc)
"""
