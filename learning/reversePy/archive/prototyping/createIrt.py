import sqlite3
import time
import statistics

import numpy as np

#Get time to calulate program time
start = time.time()

#Connect to database
conn = sqlite3.connect('db.db')

print ("Pulling Data from Database")
c = conn.cursor()
c.execute("SELECT * FROM PeptideTable")
pepTable = c.fetchall()

f = open("irt.irt.csv", "w+")
f.write("sequence,irt\n")
current = ""
clist = []

#Iterate through table and pull information about each peptide and its scans
for ind in pepTable:
#for j in range(15):

    #Match peptide in table to peptide in Spectra Table
    pepID = (str(ind[0]), )
    #pepID = (str(pepTable[j][0]), )


    peptide = ind[1].split('.')[1]
    #peptide = pepTable[j][1].split('.')[1]

    if '(' in peptide or 'Z' in peptide or 'B' in peptide or 'U' in peptide:
        continue

    rt = float(ind[7])
    #rt = pepTable[j][7]

    #f.write(peptide + "," + str(rt) +"\n")

    if peptide != current and current != "":
        rtime = statistics.median(clist)
        f.write(current + "," + str(rtime) +"\n")
        clist = []
        current = peptide
        clist.append(rt)
    else:
        current = peptide
        clist.append(rt)
f.close()

print ("Finished Pulling Data from Database")
pullData = time.time()
print ("Time for section: " + str(round(pullData - start)))
print ("Elapsed Time: " + str(round(pullData - start)) + "\n")


