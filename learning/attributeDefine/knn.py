import sqlite3
import time
import math
import numpy as np
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import NearestNeighbors
#Get time to calulate program time
start = time.time()

features = []
labels = []
data = []

base = "data/D2_pPro_B-0"
ext = ".sqt"

for i in range(1, 10):
    #file = open(inputPath, "r")
    file = open(base + str(i) + ext, "r")

    done = False
    mBool = False # Boolean for in M bounds
    normal = False

    maxIntensity = 0
    minIntensity = 10000
    maxNum = 0
    minNum = 100000
    maxLength = 0
    minLength = 1000000

    for line in file:
        line = line.split("\t")
        if line[0] == "S":
            chargeState = int(line[3])
            m1 = float(line[6])
            intensity = float(line[7])
            if intensity > maxIntensity:
                maxIntensity = intensity
            elif intensity < minIntensity:
                minIntensity = intensity
            pTrue = float(line[8])
            numMatched = int(line[9])
            if numMatched > maxNum:
                maxNum = numMatched
            elif numMatched < minNum:
                minNum = numMatched
            data.append(chargeState)
            data.append(intensity)
            data.append(pTrue)
            data.append(numMatched)
            done = False
        elif line[0] == "M" and mBool == False and done == False:
            m2 = float(line[3])
            xCorr = float(line[5])
            peptide = line[9]
            lenPeptide = len(peptide.split('.')[1])
            if lenPeptide > maxLength:
                maxLength = lenPeptide
            elif lenPeptide < minLength:
                minLength = lenPeptide
            normalized = abs((m1 - m2)/m1)
            normalized = math.modf(normalized)[0]

            #data.append(m2)
            data.append(xCorr)
            data.append(lenPeptide)
            data.append(normalized)
            mBool = True

        elif line[0] != "M" and mBool == True and done == False:
            if not (line[1][0] == "R" and line[1][1] == "e" and line[1][7] == "_"):
                normal = True
        elif line[0] == "M" and mBool == True and done == False:
            data.append(line[4])
            if normal == True:
                labels.append(1)
                #labels.append(1)
            else:
                labels.append(0)
                #labels.append(0)
            features.append(data)

            data = []
            normal = False
            mBool = False
            done = True

print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")


features = np.array(features)
labels = np.array(labels)
print(features.shape, features.dtype)

X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.01, random_state=0)
neighbors = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(X_train)
distances, indices = neighbors.kneighbors(X_test)


correct = 0
wrong = 0
cnt = 0
pred0 = 0

for i, point in enumerate(indices):

    #Label Point
    target = True
    Fcnt = 0
    Tcnt = 0
    for neighborIndex in point:
        if y_train[neighborIndex] == 0:
            Fcnt += 1
        else:
            Tcnt += 1

    if Fcnt >= 1:
        target = False
    else:
        target = True

    #Check Point
    if target:
        if y_test[i] == 1:
            correct += 1
            cnt += 1
        elif y_test[i] == 0:
            wrong += 1
            cnt += 1
    else:
        pred0 += 1
print("Correct: ", correct)
print("Wrong: ", wrong)
print("Total", cnt)

print("Inconclusive: ", pred0)


writeData = time.time()
print ("Time for section: " + str(round(writeData - pullData)))
print ("Elapsed Time: " + str(round(writeData - start)) + "\n")
