# libsearch

To build consensus library:
java -cp /path/to/classpath/ edu.scripps.dia.LibraryBuilder /path/to//DTASelect-filter.txt /path/to//search.xml  /path/to/target.ms2

DTASelect-filter.txt: DTASelect output file

search.xml - search param file

target.ms2 -ouput file 

To Run:
 java -cp /path/to/classpath/ edu.scripps.dia.ConsensusLibrarySearchEngine /path/to/source.ms2 /path/to/search.xml /path/to/target.ms2 /path/to/output.txt
 
source.ms2 - spectra file to search against

search.xml - search param file

target.ms2 - file generated from building consensus library 
