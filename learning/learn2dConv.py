import sys
import pickle
import time
import random
import math
import gc

import tensorflow as tf
from tensorflow import keras

import numpy as np
import matplotlib.pyplot as plt
start = time.time()

#PROGRAM VARS

Y_TOP_BOUND  = 100000
yBins = 10
BIN_WIDTH_Y = Y_TOP_BOUND/yBins

#BIN_WIDTH_Y = 10
BIN_WIDTH_X = 1
NUM_EPOCHS = 25
NUMBER_TO_TEST = 100

print("Running program")
print("Training Data with ", NUM_EPOCHS, " epochs.")
print("Testing accuracy on ", NUMBER_TO_TEST, " elements.")
print("\n")

'''--------------------------------------------------------------------------'''

print("Gathering Data from Files")
with open ('data/spectrum', 'rb') as sp:
    spectrums = pickle.load(sp)
with open ('data/labels', 'rb') as lp:
    labels = pickle.load(lp)
with open ('data/indexList', 'rb') as lp:
    indexList = pickle.load(lp)

#Shuffle Lists
shuffle = list(zip(spectrums, labels, indexList))
random.shuffle(shuffle);
spectrums, labels, indexList = list(zip(*shuffle));
spectrums = list(spectrums);
labels = list(labels);
indexList = list(indexList);

testSpec = []
for x in range(0, NUMBER_TO_TEST):
    testSpec.append(spectrums[x])

outputList = list(set(indexList));
outputList2 = list(set(labels));

'''--------------------------------------------------------------------------'''
#Find max x and y1
MAX_X = 0;
MAX_Y = 0;
for j in spectrums:
    for i in j:
        if i[0] > MAX_X:
            MAX_X = i[0]
        if i[1] > MAX_Y:
            MAX_Y = i[1]


MAX_X = int(math.ceil(MAX_X / 1000.0)) * 1000
MAX_Y =int(math.ceil(MAX_Y / 1000.0)) * 1000

xBins = math.ceil(MAX_X / BIN_WIDTH_X);
yBins = math.ceil(Y_TOP_BOUND/ BIN_WIDTH_Y);
#yBins = math.ceil(MAX_Y / BIN_WIDTH_Y);

totalBins = xBins + xBins * (yBins-1)
print( "Putting all coords into ", totalBins, " bins.")
print( "XBINS: ", xBins, " YBINS: ", yBins)
print( "Bins are ", BIN_WIDTH_X, " by ", BIN_WIDTH_Y)
print( "Max X: ", MAX_X, " Max Y: ", MAX_Y, "\n")
#Manipulate the Data

#Initialize a binList with 0's

binArray = []
binList = np.zeros((int(xBins), int(yBins)))
yBinTemp = []

for ind in range(0, len(spectrums)):
    print(ind,"/",len(spectrums), end="\r", flush = True)
    numPoints = len(spectrums[ind])
    for cnt, i in enumerate(spectrums[ind]):

        """print("-----------------------------------------")
        print("MAX_X: ", MAX_X, " MAX_Y: ", MAX_Y)
        print("Number Points: ", numPoints)
        print("BINS: ", xBins, " ", yBins)
        print ("BIN Width: ", BIN_WIDTH_X, BIN_WIDTH_Y)
        print("Binlength: ", len(binList))
        print("x val: ", i[0])
        print("y val: ", i[1])"""


        if i[0] == 0 and i[1] == 0:
            numPoints -= 1
        else:
            x = i[0]
            y = i[1]
            if  y > Y_TOP_BOUND:
                y = Y_TOP_BOUND - 1
            x = math.floor(x/BIN_WIDTH_X);
            y = math.floor(y/BIN_WIDTH_Y);
            #binNumber = int(x + (y*xBins));
            #print(x, " | ", y)
            #print(binNumber)
            binList[x,y] += 1
            #print(binList[x][y])
            """    for i in binList:
                    i = float(i/float(numPoints) * 50)
                    i = i/2"""
    np.divide(binList,numPoints)
    #binCopy = binList.copy()
    binArray.append(binList)
    #np.concatenate((binArray, binList))


    #binList.clear()
    np.zeros((xBins, yBins))
binArray = np.array(binArray)
binArray = np.reshape(binArray, (len(indexList), xBins, yBins, 1))
print(binArray.shape)
testSet = []
testLab = []
gc.collect()
#Create test set
for i in range(0, NUMBER_TO_TEST):
    testSet.append(binArray[i])
    testLab.append(indexList[i])
    np.delete(indexList, i)
    np.delete(binArray, i)


print ("Finished Putting into Bins")
binTime = time.time()
print ("Time for section: " + str(round(binTime - start)))
print ("Elapsed Time: " + str(round(binTime - start)) + "\n")
gc.collect()
#Create test set
for i in range(0, NUMBER_TO_TEST):
    testSet.append(binArray[i])
    testLab.append(indexList[i])
    np.delete(indexList, i)
    np.delete(binArray, i)


print ("Finished Putting into Bins")
binTime = time.time()
print ("Time for section: " + str(round(binTime - start)))
print ("Elapsed Time: " + str(round(binTime - start)) + "\n")
'''--------------------------------------------------------------------------'''

inputs = len(spectrums)
inputLayers = len(binArray)
outputLayers = len(set(indexList))


outputList = np.array(indexList)
testSet = np.array(testSet)
testLab = np.array(testLab)


print("Inputs (spectrums): ",inputs)
print("\tShapes: ", binArray.shape);
print("\tType: ", binArray.dtype);
print("\tInput Layers: ", inputLayers, "\n")

print("Outputs:", )
print("\tShapes: ", outputList.shape);
print("\tType: ", outputList.dtype);
print("\tOutput Layers: ", outputLayers)

"""print("Training Inputs: ")
print("\tShapes: ", testSet.shape);
print("\tType: ", testSet.dtype);
print("\tInput Layers: ", inputLayers, "\n")
print("Outputs:", )
print("\tShapes: ", testLab.shape);
print("\tType: ", testLab.dtype);
print("\tOutput Layers: ", outputLayers)
print("\n")"""



model = keras.Sequential([
    keras.layers.Conv2D(32, (3,3), activation=tf.nn.relu, input_shape = (xBins, yBins, 1)),
    keras.layers.Conv2D(64, (3,3), activation=tf.nn.relu),
    keras.layers.Flatten(),
    keras.layers.Dense(inputLayers/2, activation=tf.nn.relu),
    keras.layers.Dense(outputLayers, activation=tf.nn.softmax)
])

model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])


'''--------------------------------------------------------------------------'''

#Train the data
model.fit(binArray, indexList, epochs = NUM_EPOCHS)

result = model.predict(testSet)

for i, ele in enumerate(result):
    maxVal = 0;
    val2 = 0;
    val3 = 0;
    ind = 0
    ind2 = 0
    ind3 = 0
    for k, g in enumerate(ele):
        if g > maxVal:
            val3 = val2
            ind3 = ind2
            val2 = maxVal
            ind2 = ind
            maxVal = g
            ind = k

    predicted = outputList2[ind]
    actual = outputList2[testLab[i]]
    print("Test ", i)
    print(ind, "  |  ", testLab[i])
    print("Predicted: ", predicted)
    print("ActualVal: ", actual, "\n")
    if predicted != actual:
        print(ind,": ", maxVal, "|", end='')
        print(ind2,": ", val2, "|", end='')
        print(ind3,": ", val3, "|", end='')
        print ("\n")
        print(testLab[i], ": ", ele[testLab[i]])
        print("##############################################")
    else:
        print("----------------------------------------------")
"""        x1 = [x[0] for x in testSpec[ind]]
        y1 = [x[1] for x in testSpec[ind]]
        x2 = [x[0] for x in testSpec[i]]
        y2 = [x[1] for x in testSpec[i]]
        plt.style.use('ggplot')
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.scatter(x1, y1, c='r', label = 'first')
        ax1.scatter(x2, y2, c='b', label = 'second')
        plt.legend(loc='upper left');
        plt.show()"""


#Test Acccuracy
np.array(testLab)
test_loss, test_acc = model.evaluate(testSet, testLab)
print('Test accuracy:', test_acc)


print ("Finished Training")
train = time.time()
print ("Time for section: " + str(round(train - binTime)))
print ("Elapsed Time: " + str(round(train - start)) + "\n")
