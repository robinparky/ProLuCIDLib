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

with open (inputPath + 'outputLabels', 'rb') as lp:
    outputLabels = pickle.load(lp)
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


if showResults == "True" or showResults == "true" or showResults == "t" or showResults == "T":
    for i, ele in enumerate(result):
        #print(ele)
        cnt=0; #Variables for keeping track of max values
        maxVal = 0;
        val2 = 0;
        val3 = 0;
        ind = 0
        ind2 = 0
        ind3 = 0
        #Iterate through result, setting max values
        for k, g in enumerate(ele):
            if g > maxVal:
                val3 = val2
                ind3 = ind2
                val2 = maxVal
                ind2 = ind
                maxVal = g
                ind = k
            elif g > val2:
                val3 = val2
                ind3 = ind2
                val2 = g
                ind2 = k
            elif g > val3:
                val3 = g
                ind3 = k

        predicted = outputLabels[ind]
        actual = testLab[i]
        print("Test ", i)
        print(ind, "  |  ", testInd[i])
        print("Predicted: ", predicted)
        print("ActualVal: ", actual, "\n")
        print(ind,": ", maxVal, "|", end='')
        print(ind2,": ", val2, "|", end='')
        print(ind3,": ", val3, "|", end='')
        print ("\n")

        if predicted != actual or ind != testInd[i]:
        #if 1 ==1:

            ###################################################################
            x1 = [x[0] for x in testSpec[i]]
            y1 = [x[1] for x in testSpec[i]]
            index = list(labelList).index(actual)
            actualSpectrum = spectrums[index];
            x2 = [x[0] for x in actualSpectrum]
            y2 = [x[1] for x in actualSpectrum]
            plt.style.use('ggplot')
            plt.subplot(2,2,1)
            plt.title("Spectrum Comparison")
            plt.scatter(x1, y1, c='r', label = 'Predicted: ' + str(predicted) + " | " + str(testId[i]),s=1)
            plt.scatter(x2, y2, c='b', label = 'Actual: ' + str(actual)+ " | " + str(idList[index]),s=1)
            plt.legend(loc='upper left');
            #######################################################################
            plt.subplot(2,2,2)
            xValueList1 = [];
            yValueList1 = [];
            xValueList2 = [];
            yValueList2 = [];
            predictedBins = testBins[i]
            actualBins = binArray[index]
            for j in range(0, totalBins):
                if predictedBins[j] !=0:
                    xValueList1.append(j);
                    yValueList1.append(predictedBins[j])
                if actualBins[j] !=0:
                    xValueList2.append(j);
                    yValueList2.append(actualBins[j])

            plt.scatter(xValueList1, yValueList1, c='r', label = 'Predicted: ' + str(predicted)+ " | " + str(testId[i]),s=1)
            plt.scatter(xValueList2, yValueList2, c='b', label = 'Actual: ' + str(actual) + " | " + str(idList[index]),s=1)
            plt.legend(loc='upper left')
            plt.title("Bin Comparison")
            plt.xlabel("Bins")
            plt.ylabel("Percentage(Normalized)")
            #############################################################################
            x1 = [x[0] for x in testSpec[i]]
            y1 = [x[1] for x in testSpec[i]]
            index = list(labelList).index(predicted)
            actualSpectrum = spectrums[index];
            x2 = [x[0] for x in actualSpectrum]
            y2 = [x[1] for x in actualSpectrum]

            plt.subplot(2,2,3)
            plt.title("Spectrum Comparison")
            plt.scatter(x1, y1, c='r', label = 'Predicted: ' + str(predicted)+ " | " + str(testId[i]),s=1)
            plt.scatter(x2, y2, c='b', label = 'Actual: ' + str(labelList[index])+ " | " + str(idList[index]),s=1)
            plt.legend(loc='upper left');
            #######################################################################
            plt.subplot(2,2,4)
            xValueList1 = [];
            yValueList1 = [];
            xValueList2 = [];
            yValueList2 = [];
            predictedBins = testBins[i]
            actualBins = binArray[index]
            for j in range(0, totalBins):
                if predictedBins[j] !=0:
                    xValueList1.append(j);
                    yValueList1.append(predictedBins[j])
                if actualBins[j] !=0:
                    xValueList2.append(j);
                    yValueList2.append(actualBins[j])
            plt.scatter(xValueList1, yValueList1, c='r', label = 'Predicted: ' + str(predicted)+ " | " + str(testId[i]),s=1)
            plt.scatter(xValueList2, yValueList2, c='b', label = 'Actual: ' + str(labelList[index])+ " | " + str(idList[index]),s=1)
            plt.legend(loc='upper left')
            plt.title("Bin Comparison")
            plt.xlabel("Bins")
            plt.ylabel("Percentage(Normalized)")

            """manager = plt.get_current_fig_manager()
            manager.resize(*manager.window.maxsize())
            path = str("falsePredict/"+ str(cnt) + ".png")
            cnt +=1;
            plt.savefig(path)"""
            print("##############################################")
            plt.show()

        else:
            print("----------------------------------------------")

#Test Acccuracy
np.array(testInd)
test_loss, test_acc = model.evaluate(testBins, testInd)
print('Test accuracy:', test_acc)
