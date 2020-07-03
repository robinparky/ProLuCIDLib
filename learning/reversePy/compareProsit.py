import pandas as pd
import copy
import sys
import json

if len(sys.argv) != 4:
    print("python compareProsit.py prositPrediction.csv ourPrediction.json experimental.json"

prositDf  = pd.read_csv(sys.argv[1])

templateObject = {'y1':[], 'b1':[], 'y2':[], 'b2':[]}


currentPeptide = ""

prositDict = {}

tmpObj = {'y1':[], 'b1':[], 'y2':[], 'b2':[]}
for index, row in prositDf.iterrows():
    #Check if new peptide
    if index == 0:
        currentPeptide = row["StrippedPeptide"]
    elif(row["StrippedPeptide"] != currentPeptide):

        #print("-----------")
        #print(tmpObj['y1'])
        #print(tmpObj['y2'])
        #print(tmpObj['b1'])
        #print(tmpObj['b2'])
        prositDict[currentPeptide] = copy.copy(tmpObj)
        currentPeptide = row["StrippedPeptide"]
        tmpObj = {'y1':[], 'b1':[], 'y2':[], 'b2':[]}

    #print(currentPeptide)

    fragType = str(row["FragmentType"])
    fragCharge = str(row["FragmentCharge"])
    intensity = row["RelativeIntensity"]

    tmpObj[fragType+fragCharge].append(intensity)

"""
for key, value in prositDict.items():
    print(key)
    print(value['y1'])
    print(value['b1'])
    print(value['y2'])
    print(value['b2'])
    sys.exit()
"""
# JSON file

predictions = open (sys.argv[2], "r")
experimental = open(sys.argv[3], "r")

# Reading from file
predData = json.loads(predictions.read())
expData = json.loads(experimental.read())


for key, value in prositDict.items():
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    for fragtype, arr in value.items():
        print(fragtype, key)
        exp = expData[key][fragtype]
        pred = predData[key][fragtype]
        prosit = value[fragtype]

        while len(prosit) != len(exp):
            prosit.append(999999)

        #print(len(exp))
        #print(len(pred))
        #print(len(prosit))

        print('{:15s} {:15s} {:15s}'.format("Experimental","Prediction","Prosit"))
        for i in range(len(exp)):
            print('{:15f} {:15f} {:15f}'.format(float(exp[i]),float(pred[i]),float(prosit[i])))
