import sys
import pickle
import time
import random

import tensorflow as tf
from tensorflow import keras

import numpy as np
import matplotlib.pyplot as plt

start = time.time()
print("Gathering Data from Files")
with open ('data/spectrumsList', 'rb') as sp:
    spectrums = pickle.load(sp)
with open ('data/labelList', 'rb') as lp:
    labels = pickle.load(lp)
with open ('data/indexList', 'rb') as lp:
    indexList = pickle.load(lp)



MAX_X = 0;
MAX_Y = 0;

for j in spectrums:
    for i in j:
        if i[0] > MAX_X:
            MAX_X = i[0]
        if i[1] > MAX_Y:
            MAX_Y = i[1]



numInputs = len(spectrums)
numOutputs = len(set(indexList))

print("\tSpectrums: ", len(spectrums))
print("\tElement Data Points",len(spectrums[0]))
print("\tUnique Labels: ", numOutputs)

spectrums[0] = sorted(spectrums[0])

#Sort all spectrums in ascending order.
for i, ele in enumerate(spectrums):
    spectrums[i] = sorted(spectrums[i])
plt.style.use('ggplot')
#PRINT ALL DIFFERENT LABELS
var = input("Show quantity Labels? ")
if var == "y":
    plt.title("Quantity of Labels")
    plt.hist(labels, len(labels))
    plt.ylabel("Count")
    plt.xlabel("label", )
    plt.xticks(labels, rotation='vertical')
    plt.show()



#PRINT ALL UNIQUE Label and all corresponding points
var = input("Show all points for each unique Label? ")
if var == "y":
    for cnt, i in enumerate(set(labels)):
        indlist = [y for y, x in enumerate(labels) if x == i];
        speclist = []
        for k in indlist:
            speclist.append(spectrums[k]);
        x = []
        y = []
        for j in speclist:
            xi = [x[0] for x in j]
            yi = [x[1] for x in j]
            x = x + xi
            y = y + yi
        plt.title(str(i) + " with " + str(len(indlist)) + " Instances")
        plt.scatter(x, y)
        plt.xlim(0, MAX_X  + (MAX_X * .1))
        plt.ylim(0, MAX_Y + (MAX_Y * .1))
        plt.show()

var = input("Iterate through all spectrums 1 label at a time? ")
if var == "y":
    for cnt, i in enumerate(set(labels)):
        indlist = [y for y, x in enumerate(labels) if x == i];
        speclist = []
        for k in indlist:
            speclist.append(spectrums[k]);
        for j in speclist:
            x = [x[0] for x in j]
            y = [x[1] for x in j]
            plt.title(str(i) + " with " + str(len(j)) + " Seperate Points")
            plt.scatter(x, y)
            plt.xlim(0, MAX_X  + (MAX_X * .1))
            plt.ylim(0, MAX_Y + (MAX_Y * .1))
            plt.show()

#PRINT ALL DATA POINTS ON 1 GRAPH
var = input("Show all data points on 1 graph? ")
if var == "y":
    x = []
    y = []
    for i in spectrums:
        xi = [x[0] for x in i]
        yi = [x[1] for x in i]
        x = x + xi
        y = y + yi
        print(xi, ", ", yi)
    plt.scatter(x,y)
    plt.show()
var = input("Search for spectrum:");
individual = input("group or ind?")
if individual == "group":
    for cnt, i in enumerate(set(labels)):
        if i == var:
            indlist = [y for y, x in enumerate(labels) if x == i];
            speclist = []
            for k in indlist:
                speclist.append(spectrums[k]);
            x = []
            y = []
            for j in speclist:
                xi = [x[0] for x in j]
                yi = [x[1] for x in j]
                x = x + xi
                y = y + yi
            plt.title(str(i) + " with " + str(len(indlist)) + " Instances")
            plt.scatter(x, y)
            plt.xlim(0, MAX_X  + (MAX_X * .1))
            plt.ylim(0, MAX_Y + (MAX_Y * .1))
            plt.show()
elif individual == "ind":
    for cnt, i in enumerate(set(labels)):
        if i == var:
            indlist = [y for y, x in enumerate(labels) if x == i];
            speclist = []
            for k in indlist:
                speclist.append(spectrums[k]);
            for j in speclist:
                x = [x[0] for x in j]
                y = [x[1] for x in j]
                plt.title(str(i) + " with " + str(len(j)) + " Seperate Points")
                plt.scatter(x, y)
                plt.xlim(0, MAX_X  + (MAX_X * .1))
                plt.ylim(0, MAX_Y + (MAX_Y * .1))
                plt.show()
else:
    print("Invalid response")
