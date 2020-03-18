import sqlite3
import sys
from pyteomics import mass as massC
import struct
import numpy as np
import pandas as pd

import time

#t
#Connect to database
conn = sqlite3.connect('/home/bernard/full.db')
#conn = sqlite3.connect('~/full.db')
#conn = sqlite3.connect('projects2012_01_05_20_21166.db')
#conn = sqlite3.connect('testLibDuplicateSpectra.db')
#conn = sqlite3.connect(DATABASE_PATH)

print ("Pulling Data from Database")
t0 = time.time()
c = conn.cursor()
c.execute("SELECT * FROM PeptideTable")
pepTable = c.fetchall()

cnt = 0

#Iterate through table and pull information about each peptide and its scans
for ind in pepTable:

    pepID = (str(ind[0]), )
    c.execute('SELECT *,rowid FROM SpectraTable WHERE peptideID=?', pepID)
    spectrums = c.fetchall()
t1 = time.time()
print("Time taken: ", t1 - t0)
