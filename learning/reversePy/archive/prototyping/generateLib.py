import json
import sqlite3
from pyteomics import mass

c2 = ""
c3 = ""
irt = ""

with open('pepList_charge2.prediction.ions.json') as jsonFile:
    c2 = json.load(jsonFile)

with open('pepList_charge3.prediction.ions.json') as jsonFile:
    c3 = json.load(jsonFile)


conn = sqlite3.connect('library.db')
c = conn.cursor()
c.execute("CREATE TABLE IF NOT EXISTS predictions(Sequence varchar,\
                                                  PrecursorMZ float,\
                                                  Charge integer,\
                                                  Retention_Time float,\
                                                  B_Frag_Index varchar,\
                                                  B_Frag_Int varchar,\
                                                  B_Frag_Loss varchar,\
                                                  Y_Frag_Index varchar,\
                                                  Y_Frag_Int varchar,\
                                                  Y_Frag_Loss varchar)")
peptide = ""
rt = 0

f = open('pepList.prediction.irt.csv', 'r')
for ele in zip(f, c2, c3):
    line = ele[0].split(',')
    peptide = line[0]
    rt = line[1]

    preMz2 = mass.calculate_mass(sequence = peptide, ion_type = 'b', charge=2)
    c2 = ele[1]

    pep2 = ele[1]['peptide']
    pep3 = ele[2]['peptide']
    if peptide != pep2 or peptide != pep3 or pep2 != pep3:
        continue

    bidx = ''
    bint = ''
    bloss = ''
    yidx = ''
    yint = ''
    yloss = ''

    for key, mzarr in c2['ions'].items():
        for mz in mzarr:
            if mz != 0:
                if key == 'b1':
                    bloss+=('Null,')
                    bidx+=(str(mzarr.index(mz) + 1) + ",")
                    bint+=(str(mz)+ ",")
                elif key == 'y1':
                    yloss+=('Null,')
                    yidx+=(str(mzarr.index(mz) + 1) + ",")
                    yint+=(str(mz)+ ",")
                elif key == 'bn1':
                    bloss+=('nh3,')
                    bidx+=(str(mzarr.index(mz) + 1) + ",")
                    bint+=(str(mz)+ ",")
                elif key == 'yn1':
                    yloss+=('nh3,')
                    yidx+=(str(mzarr.index(mz) + 1) + ",")
                    yint+=(str(mz)+ ",")
                elif key == 'bo1':
                    bloss+=('h2o,')
                    bidx+=(str(mzarr.index(mz) + 1) + ",")
                    bint+=(str(mz)+ ",")
                elif key == 'yo1':
                    yloss+=('h2o,')
                    yidx+=(str(mzarr.index(mz) + 1) + ",")
                    yint+=(str(mz)+ ",")
    """
    print(peptide)
    print(rt)
    print(preMz2)
    print(bidx)
    print(bint)
    print(bloss)
    print(yidx)
    print(yint)
    print(yloss)
    break
    """
    c.execute("INSERT INTO predictions(Sequence, PrecursorMZ, Charge, Retention_Time,\
                                       B_Frag_Index, B_Frag_Int, B_Frag_Loss,\
                                       Y_Frag_Index, Y_Frag_Int, Y_Frag_Loss)\
               VALUES(?,?,?,?,?,?,?,?,?,?)", (peptide, preMz2, 2, rt, bidx,bint,bloss,yidx,yint,yloss))
    conn.commit()

    preMz3 = mass.calculate_mass(sequence = peptide, ion_type = 'b', charge=3)
    c3 = ele[2]

    bidx = ''
    bint = ''
    bloss = ''
    yidx = ''
    yint = ''
    yloss = ''

    for key, mzarr in c3['ions'].items():
        for mz in mzarr:
            if mz != 0:
                if key == 'b1':
                    bloss+=('Null,')
                    bidx+=(str(mzarr.index(mz) + 1) + ",")
                    bint+=(str(mz)+ ",")
                elif key == 'y1':
                    yloss+=('Null,')
                    yidx+=(str(mzarr.index(mz) + 1) + ",")
                    yint+=(str(mz)+ ",")
                elif key == 'bn1':
                    bloss+=('nh3,')
                    bidx+=(str(mzarr.index(mz) + 1) + ",")
                    bint+=(str(mz)+ ",")
                elif key == 'yn1':
                    yloss+=('nh3,')
                    yidx+=(str(mzarr.index(mz) + 1) + ",")
                    yint+=(str(mz)+ ",")
                elif key == 'bo1':
                    bloss+=('h2o,')
                    bidx+=(str(mzarr.index(mz) + 1) + ",")
                    bint+=(str(mz)+ ",")
                elif key == 'yo1':
                    yloss+=('h2o,')
                    yidx+=(str(mzarr.index(mz) + 1) + ",")
                    yint+=(str(mz)+ ",")
    c.execute("INSERT INTO predictions(Sequence, PrecursorMZ, Charge, Retention_Time,\
                                       B_Frag_Index, B_Frag_Int, B_Frag_Loss,\
                                       Y_Frag_Index, Y_Frag_Int, Y_Frag_Loss)\
               VALUES(?,?,?,?,?,?,?,?,?,?)", (peptide, preMz3, 3, rt, bidx,bint,bloss,yidx,yint,yloss))
    conn.commit()
print("Exiting Successfully")
