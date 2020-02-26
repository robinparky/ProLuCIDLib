import sqlite3
import sys
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd
import time

if len(sys.argv) != 2:
    print("Error with command line inputs")
    sys.exit(0)
else:
    INPUT_PATH = sys.argv[1]


#Encode Sequence
df = pd.read_pickle(INPUT_PATH)
#read the matrix a csv file on github

outFile = open("fullTest.peptide.csv", "w")
outFile.write("sequence, fullSequence, scan, fileName\n")
for i in range(100):
    target = df.iloc[i]
    outFile.write(target["peptide"]+"," + target['peptideFull'] + "," + str(target['scan'])+ "," + target['fileName'] +  "\n")
outFile.close()
