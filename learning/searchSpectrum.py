import sqlite3
import matplotlib.pyplot as plt
import struct

def convertFloat(element):
    return struct.iter_unpack('>f', element)

peptide = input("Enter peptide to search for: ")
ind = input("Print 2 at Spectrums at a time?")
allSpecs = [];
oldSpec = [];

init = 0

databasePath = 'data/testLibDuplicateSpectra.db'
conn = sqlite3.connect(databasePath)

c = conn.cursor()
c.execute("SELECT * FROM PeptideTable WHERE sequenceCS=?", (peptide, ))
peptideList = c.fetchone()
c.execute("SELECT * FROM SpectraTable WHERE peptideID=?", (str(peptideList)[1], ))
spectrums = c.fetchall();
if len(spectrums) != 0:
    for element in spectrums: #Search returned list of matched spectrums
        spectrum = [];
        mzArr = convertFloat(element[1])
        intArr = convertFloat(element[2])
        for i, j in zip(mzArr, intArr):
            coords = []
            mz = float(i[0])
            intensity = float(j[0]);
            coords.append(mz);
            coords.append(intensity);
            spectrum.append(coords);
        allSpecs.append(spectrum)
        if ind == "y":
            x1 = [];
            y1 = [];
            x2 = [];
            y2 = [];
            plt.style.use('ggplot')
            plt.title("Spectrum Comparison")
            x1 = [x[0] for x in spectrum]
            y1 = [x[1] for x in spectrum]
            x2 = [x[0] for x in oldSpec]
            y2 = [x[1] for x in oldSpec]
            plt.scatter(x1, y1, c='r',s=1)
            plt.scatter(x2, y2, c='b',s=1)
            plt.show()
            oldSpec = spectrum
    x = []
    y = []
    for j in allSpecs:
        xi = [x[0] for x in j]
        yi = [x[1] for x in j]
        x = x + xi
        y = y + yi
    plt.title(peptide + " with " + str(len(allSpecs)) + " Instances")
    plt.scatter(x, y)
    plt.show()
