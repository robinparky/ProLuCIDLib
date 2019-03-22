import sqlite3
import matplotlib.pyplot as plt
import struct

def convertFloat(element):
    return struct.iter_unpack('>f', element)

id1 = input("Enter SpectrumID 1: ")
id2 = input("Enter SpectrumID 2: ")
spectrum1 = [];
spectrum2 = [];

databasePath = 'data/testLibDuplicateSpectra.db'
conn = sqlite3.connect(databasePath)
c = conn.cursor()

c.execute("SELECT * FROM SpectraTable WHERE rowid=?", (id1, ))
spectrums = c.fetchone();
spectrum = [];
mzArr = convertFloat(spectrums[1])
intArr = convertFloat(spectrums[2])
for i, j in zip(mzArr, intArr):
    coords = []
    mz = float(i[0])
    intensity = float(j[0]);
    coords.append(mz);
    coords.append(intensity);
    spectrum.append(coords);
spectrum1 = spectrum

c.execute("SELECT * FROM SpectraTable WHERE rowid=?", (id2, ))
spectrums = c.fetchone();
spectrum = [];

mzArr = convertFloat(spectrums[1])
intArr = convertFloat(spectrums[2])
for i, j in zip(mzArr, intArr):
    coords = []
    mz = float(i[0])
    intensity = float(j[0]);
    coords.append(mz);
    coords.append(intensity);
    spectrum.append(coords);
spectrum2 = spectrum


x1 = [];
y1 = [];
x2 = [];
y2 = [];
plt.style.use('ggplot')
plt.title("Spectrum Comparison")
x1 = [x[0] for x in spectrum1]
y1 = [x[1] for x in spectrum1]
x2 = [x[0] for x in spectrum2]
y2 = [x[1] for x in spectrum2]
plt.scatter(x1, y1, c='r',s=1)
plt.scatter(x2, y2, c='b',s=1)
plt.show()
oldSpec = spectrum
