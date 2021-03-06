from __future__ import division
import sqlite3
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

if len(sys.argv) != 5:
    print("Error with command line inputs")
    sys.exit(0)
else:
    inputPath = sys.argv[1]
    neuralPath = sys.argv[2]
    showResults = sys.argv[3]
    databasePath = sys.argv[4]

'''--------------------------------------------------------------------------'''
print("Gathering Data from Files")
with open (inputPath + 'spectrumsList', 'rb') as sp:
    spectrums = pickle.load(sp)
sp.close()
with open (inputPath + 'labelList', 'rb') as lp:
    labelList = pickle.load(lp)
sp.close()
with open (inputPath + 'indexList', 'rb') as lp:
    indexList = pickle.load(lp)
sp.close()
with open (inputPath + 'idList', 'rb') as lp:
    idList = pickle.load(lp)
sp.close()
with open (inputPath + 'binArray', 'rb') as lp:
    binArray = pickle.load(lp)
sp.close()
with open (inputPath + 'testSpec', 'rb') as sp:
    testSpec = pickle.load(sp)
sp.close()
with open (inputPath + 'testBins', 'rb') as lp:
    testBins = pickle.load(lp)
sp.close()
with open (inputPath + 'testInd', 'rb') as lp:
    testInd = pickle.load(lp)
sp.close()
with open (inputPath + 'testLab', 'rb') as lp:
    testLab = pickle.load(lp)
sp.close()
with open (inputPath + 'testId', 'rb') as lp:
    testId = pickle.load(lp)
sp.close()

with open (inputPath + 'testMass', 'rb') as lp:
    testMass = pickle.load(lp)
sp.close()

with open (inputPath + 'testRt', 'rb') as lp:
    testRt = pickle.load(lp)
sp.close()

with open (inputPath + 'outputLabels', 'rb') as lp:
    outputLabels = pickle.load(lp)
sp.close()
'''--------------------------------------------------------------------------'''

inputs = len(spectrums)
totalBins = len(binArray[0])
inputLayers = len(binArray)
output = np.unique(indexList)
outputLayers = len(outputLabels)

print("inputs", inputs)
print("outputs", outputLayers)


def create_model():
    model = keras.Sequential()
    with tf.device('/gpu:0'):
        model.add(keras.layers.InputLayer(input_shape = (totalBins, )))
        model.add(keras.layers.Dense(1000, activation=tf.nn.relu))
    with tf.device('/gpu:1'):
        model.add(keras.layers.Dense(1000, activation=tf.nn.relu))
        model.add(keras.layers.Dense(1000, activation=tf.nn.relu))
    with tf.device('/gpu:2'):
        model.add(keras.layers.Dense(1000, activation=tf.nn.relu))
        model.add(keras.layers.Dense(1000, activation=tf.nn.relu))
    with tf.device('/gpu:3'):
        model.add(keras.layers.Dense(1000, activation=tf.nn.relu))
        model.add(keras.layers.Dense(outputLayers, activation=tf.nn.softmax))
    return model

model = create_model()

model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])

model.load_weights(neuralPath)
#model = tf.keras.models.load_model(neuralPath)

result = model.predict(testBins)



conn = sqlite3.connect(databasePath)
c = conn.cursor()


indexes = []
for i in range(len(outputLabels)):
    indexes.append(i)


correctCnt = 0
cnt = 0

for i, ele in enumerate(result):
    if(len(testSpec[i]) < 50):
        continue
    predicted = ""
    precursorMass = 0
    mass = testMass[i]
    rt = testRt[i]
    found = False

    predictList = zip(ele, outputLabels)
    predictList = sorted(predictList, key = lambda t: t[0], reverse=True)

    _, predictLabels =  list(zip(*predictList))

    for label in predictLabels:
        c.execute('SELECT * FROM PeptideTable WHERE sequenceCS = ?', (label,))
        peptide = c.fetchone()
        precursorMass =  peptide[3]
        predRt = peptide[7]

        deltaRt = abs(rt - predRt)
        deltaMass =  abs(mass - precursorMass)
        deltaDec = deltaMass % 1

        """
        if deltaRt < 5:
            if deltaMass < 4:
                for j in range(4):
                    print(" ----",abs(float(j)-deltaDec))
                    if abs(j - deltaDec) < .05:
                        print("#",deltaMass,deltaRt, peptide[14])
                        print("found")
                        predicted = label
                        found = True
                        break
                if found:
                    break
        """

        if deltaMass < 4:
            for j in range(4):
                print(" ----",abs(float(j)-deltaDec))
                if abs(j - deltaDec) < .05:
                    if deltaRt < 5:
                        print("#",deltaMass,deltaRt, peptide[14])
                        print("found")
                        predicted = label
                        found = True
                        break
            if found:
                break

    if not found:
        predicted = predictList[0][1]


    actual = testLab[i]
    actInd = list(outputLabels).index(actual)
    predInd = list(outputLabels).index(predicted)
    print("Test ", i + 1)
    print("Predicted: ", predicted)
    print("ActualVal: ", actual)
    print("Spectrum Mass: ", mass)
    print("Actual PrecursorMass: ", precursorMass)
    print("Spectrum RT: ", rt)
    print("Actual Retention: ", predRt)
    print("Predicted Spectrum ID ",testId[i])
    print("Actual Spectrum ID", idList[actInd])
    print("Sampel Predicted ID", idList[predInd],"\n")

    #Ignore Charge State
    predicted = predicted[:-1]
    actual = actual[:-1]


    if predicted == actual :
        correctCnt += 1
        print("---------------------------------------------")
    else:
        #input()
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")

    cnt += 1

print("Accuracy: ", str(correctCnt/cnt))
