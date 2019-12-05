# libsearch

Current version 0.9b
## To build  library(sqlite):

java -cp ProluCIDLib.jar edu.scripps.dia.LibraryIndexer db_name.db /path_1/DTASelect-filter.txt /path_2/DTASelect-filter.txt


DTASelect-filter.txt: DTASelect output file

search.xml - search param file

test.db -output db; If it already exists will add peptides from new dtaselect-filter.txt
## To Run Search(sqlite):

java -cp /path/to/classpath/ edu.scripps.dia.LibrarySearch /path/to/test.ms2 /path/to//search.xml \[num threads]

search.xml - database_name must contain path to sqlite path

## To build consensus library:

java -cp /path/to/classpath/ edu.scripps.dia.LibraryBuilder /path/to//DTASelect-filter.txt /path/to//search.xml  /path/to/target.ms2

DTASelect-filter.txt: DTASelect output file

search.xml - search param file

target.ms2 -ouput file 

## Split Ms2 file:

java -cp /path/to/classpath/ script.Ms2FileSplit /path/to/dir numSplit

/path/to/dir- path to directory of library files; should be the output from LibraryBuilder

numSplit - number of files the library will be spit into; 500 is recommended

## To Run:

java -cp /path/to/classpath/ edu.scripps.dia.ConsensusLibrarySearchEngine /path/to/source.ms2 /path/to/search.xml /path/to/lib  /path/to/output.txt
 
source.ms2 - spectra file to search against

search.xml - search param file

/path/to/lib - path to library directory generated from ms2split

target.ms2 - file generated from building consensus library 

Test data can be download from 

http://massive.ucsd.edu
Massive data id: MSV000081719
