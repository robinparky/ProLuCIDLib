#nohup python createDf.py /home/bernard/full.db Full/ &> programOutput/createDf.txt
nohup python createDf.py /home/bsuwirjo/full.db Full/ &> programOutput/createDf.txt

#Charge3############################
nohup python createSets.py Full/c2DF.pkl Full/charge2/ 1 &> programOutput/createSetsC2.txt
#nohup python trainMs3.py Full/charge2/ 1024 500 &> programOutput/trainC2Ms2.txt
nohup python trainRt.py Full/charge2/ 1024 500 &> programOutput/trainC2Rt.txt

#nohup python predictFromFasta.py /home/bernard/peptideList.txt fullDb.db &> programOutput/predictC2.txt
#nohup python predictFromFasta.py /home/bsuwirjo/peptideList.txt fullDb.db &> programOutput/predictC2.txt

#Charge3#############################
#nohup python createSets.py Full/c3DF.pkl Full/charge3/ 1 &> programOutput/createSetsC3.txt
#nohup python trainMs3.py Full/charge3/ 1024 500 &> programOutput/trainC3Ms2.txt
#nohup python trainRt.py Full/charge3/ 1024 500 &> programOutput/trainC3Rt.txt

#nohup python predictFromFasta.py /home/bernard/peptideList.txt fullDb.db &> programOutput/predictC3.txt
#nohup python predictFromFasta.py /home/bsuwirjo/peptideList.txt fullDb.db &> programOutput/predictC3.txt
