/**
 * @file SpectrumReader.java
 * This is the source file for edu.scripps.pms.util.spectrum.SpectrumReader
 * @author Tao Xu
 * @date $Date: 2010/10/29 04:16:16 $
 */



package edu.scripps.dia.util;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

//import java.util.ListIterator;

public class SpectrumReader {

    public static final char FIRSTCHAROFHLINE = 'H';

    // delimiter of m2z and intensity
    //public static final String MZINTENSITYDELIMITER = " ";
    public static final String MZINTENSITYDELIMITER = "\\s+";
    private BufferedReader br; // reads from the ms file
    // the format of ms data file, e.g., ms1, ms2, dta
    private String fileFormat;
    private String msFileName;

    // if the spectraDelimiter is "S", then isNewFormat is true,
    // if the spectraDelimiter is ":", then isNewFormat is false
    private boolean isNewFormat = true;
    // delimiter for spectra, could be ':' or 'S', depends on the ms file type
    private char spectraDelimiter = 'S';
    
    private List<String> hlines = null;
    private String lastLine = ""; // remember the last line read
    private int numOfSpectra=0;
    
    private double [] retentiontimes; // assume no more than 100000 spectrum per file
    private double [] ioninjectiontimes; // assume no more than 100000 spectrum per file
    public SpectrumReader(String msFileName, String fileFormat) throws IOException {
        
        br = new BufferedReader(new FileReader(msFileName), 4096);
        this.fileFormat = fileFormat;
        this.msFileName = msFileName;
        readHlines();
        if (lastLine != null) {
            spectraDelimiter = lastLine.charAt(0);
//System.out.println("spectraDelimiter: " + spectraDelimiter);
            if (spectraDelimiter == 'S') {
                isNewFormat = true;
            } else if (spectraDelimiter == ':') {
                isNewFormat = false;
            } else {
                throw new FileFormatUnknownException();
            }
        }
    }
    public double [] getRetentionTimes() throws IOException {
        if(retentiontimes == null) {
            retrieveRetentionTimes();
        }
        return retentiontimes;
    }
    public double [] getIonInjectionTimes() throws IOException {
        if(ioninjectiontimes == null) {
            retrieveRetentionTimes(); // read both ion injection times and retention times
        }
        return ioninjectiontimes;
    }
    public double scanNum2RetentionTime(int scannum) throws IOException {
        return getRetentionTimes()[scannum]; 
    }
    public double scanNum2IonInjectionTime(int scannum) throws IOException {
        return getIonInjectionTimes()[scannum]; 
    }
    // retrieve both retention time and ion injection time
    private void retrieveRetentionTimes() throws IOException {
        retentiontimes = new double [100000]; // assume no more than 100000 spectrum per file
        ioninjectiontimes = new double [100000];  
        SpectrumReader sr = new SpectrumReader(msFileName, fileFormat);
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spec = it.next();
            int scannum = spec.getLoscan();
            double rettime = spec.getRetentionTime();
            double ioninjtime = spec.getIonInjectionTime();
            retentiontimes[scannum] = rettime;
            ioninjectiontimes[scannum] = ioninjtime;
        }
        sr.closeDataFile();
    }    
    public int getNumSpectra() throws IOException {
	int numSpectra = 0;
	BufferedReader br1 = new BufferedReader(new FileReader(msFileName), 4096);
	String line = null;
	while((line = br1.readLine()) != null) {
	    if(line.charAt(0) == spectraDelimiter) {
		numSpectra++;
	    }
	}
	br1.close();
	return numSpectra;
    }
    
    // this function should be called before calling getSpectra() more than once 
    public void refreshInputStream() throws IOException {
        closeDataFile();
        br = new BufferedReader(new FileReader(msFileName), 4096);
        
        readHlines(); 
        
    }
    public static void main(String args[]) throws Exception {


        String file = "/data/2/rpark/ip2_data/rpark/Jolene_Hela_loading/20141101_HeLa_1ug_BEH60_140_35_IC_DE5_5e3_1_2014_11_12_16_2015_03_19_17_30916/search/projects2016_04_01_16_95760/20141101_HeLa_1ug_BEH60_140min_35ms_IC_DE5_5e3_1.ms2";
        SpectrumReader sr = new SpectrumReader(file, "ms2");
        ArrayList<PeakList> peaklists = sr.getSpectraList();
        sr.closeDataFile();


        if(true) return;

        String line = null;
         sr = new SpectrumReader(args[0], args[1]);
	Hline h = new Hline(sr.getHlines());

        
        Iterator<PeakList> it = sr.getSpectra();
        int counter = 0;
        int numPeaks = 0;
        //boolean sortByIntensity = true;
        while (it.hasNext()) {
            PeakList list = it.next();
            //System.out.println("lowscan: " + list.getLoscan() + "\thiscan: " + list.getHiscan() + "\tnum zlines: " + list.getNumZlines());
            //System.out.println("RetentionTime: " + list.getRetentionTime());
            int scanNum = list.getLoscan();
            //System.out.println(scanNum + "\t" + sr.scanNum2RetentionTime(scanNum) + "\t" + sr.scanNum2IonInjectionTime(scanNum) + "\t" + list.numPeaks());
          //\t" + list.getHiscan() + "\t");

//	    System.out.println(list.getSpectrumWithoutHlines());
	    Peak p = null;
	    StringBuffer sb = new StringBuffer();
	//System.out.println("===>>" + list.getZlines());

	    for(Iterator<Peak> itr = list.getPeaks(); itr.hasNext();) {
		    p = itr.next();
		    sb.append(p.getM2z());
		    sb.append("\t");
		    sb.append(p.getIntensity());
		    sb.append("\n");
	    }



            for(Iterator itr = list.getZlines(); itr.hasNext(); ) {
                Zline zline = (Zline)itr.next();

               // System.out.println(zline.getM2z() + " " + zline.getChargeState());
            }


 
            //list.sort(sortByIntensity);
//            System.out.println("MaxIntens: " + list.getMaxIntensity());
//            System.out.println("MinIntens: " + list.getMinIntensity());
            //System.out.print("Min m2z: " + list.getMinM2z());
            //System.out.print("\tMax m2z: " + list.getMaxM2z());
            //System.out.print("\tRange: " + (list.getMaxM2z()-list.getMinM2z()));

            //System.out.println("\tNum peaks: " + list.numPeaks());
            numPeaks += list.numPeaks();
            counter++;
            //if (counter > 20)
            //break;
        }
        System.out.println("Total number of spectra processed: " + counter);
        System.out.println("Average number of peaks per list: " + numPeaks/counter);
        sr.closeDataFile();

        System.out.println("Finished");
    }

    public Iterator<String> getHlines() {
        return hlines.iterator();
    }
    public ArrayList<PeakList> getSpectraList() throws IOException {
        refreshInputStream();
        Iterator<PeakList> it = getSpectra();
        ArrayList<PeakList> spectraList = new ArrayList(40000);
        while(it.hasNext()) {
            spectraList.add(it.next());
        }
        return spectraList;
    }
    public Iterator<PeakList> getSpectra() throws IOException {
        return new Iterator() {
            public boolean hasNext() {
                return lastLine != null;
            }

            public Object next() {
                PeakList peaks = null;
                try {
                    peaks = getPeakList();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return peaks;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }
        };
    }

    public void closeDataFile() throws IOException {
        br.close();
    }

    private void readHlines() throws IOException {

        boolean isFirstRead = true; // if file is read for the first time
        if (hlines == null) {
            hlines = new ArrayList<String>();
        } else {
            isFirstRead = false;
        }
        while ((lastLine = br.readLine()) != null) {
            if (lastLine.charAt(0) == FIRSTCHAROFHLINE) {
                if (isFirstRead) {
                    hlines.add(lastLine);
                }
            } else if (lastLine.charAt(0) == 'S' || lastLine.charAt(0) == ':') {
                return;
            } else {
                continue;            
            }
        }
    }

    private PeakList getPeakList() throws IOException {
        PeakList peakList = new PeakList();
        peakList.setHlines(hlines);
        getListSpecificInfo(peakList);
        getPeaks(peakList);
//System.out.println("lowscan: " + peakList.getLoscan() + "\thiscan: " + peakList.getHiscan() + "\tnum zlines: " + peakList.getNumZlines());
        return peakList;
    }

    private void getPeaks(PeakList list) throws IOException {
        String[] content = null;
        while(lastLine != null && lastLine.charAt(0) != spectraDelimiter) {
	/*try {
            // process one line for m2z and intensity
            content = lastLine.split(MZINTENSITYDELIMITER);
            // add a new Peak to the peaklist
	    if (content.length < 2) {
            	list.addPeak(new Peak(Float.parseFloat(content[0]),
                                   Float.parseFloat(content[1])));
	    } else {
		list.addPeak(new Peak(Float.parseFloat(content[0]),
                                   Float.parseFloat(content[1]), Integer.parseInt(content[2])));
	    }
	} catch (Exception e) {System.out.println(e.toString() + lastLine);}
            // check next line*/
            if('H' != lastLine.charAt(0)) { 
                try {
                    // process one line for m2z and intensity
                    content = lastLine.split(MZINTENSITYDELIMITER);
                    // add a new Peak to the peaklist
                    //float m2z = Float.parseFloat(content[0]);
                    //float intens = Float.parseFloat(content[1]);
                    double m2z = Double.parseDouble(content[0]);
                    double intens = Double.parseDouble(content[1]);
                    if(m2z != 0 || intens != 0) { // avoid to add peaks that have m2z or intensity as 0
                        list.addPeak(new Peak(m2z, intens));
                    }
                } catch (Exception e) {
                    System.out.println(e.toString() + lastLine);
                }
            } else {
                // ignore H lines in the middle of the file, 
                // in case mixed instrument type was used by mistake
            }
                // check next line
            lastLine = br.readLine();
        }
    }
    private void getListSpecificInfo(PeakList list) throws IOException {

        String[] content = null;
        Zline z = null;
        if (isNewFormat){
            // parse S line
	
            content = lastLine.split("\t");
            list.setLoscan(Integer.parseInt(content[1]));
            list.setHiscan(Integer.parseInt(content[2]));
            if (fileFormat.equals("ms2")) {
                list.setPrecursorMass(Float.parseFloat(content[3]));
            }
            // parse Z lines
            while ((lastLine = br.readLine()) != null) {
                if (lastLine.charAt(0) == 'Z') {
                    content = lastLine.split("\t");
                    z = new Zline(Integer.parseInt(content[1]),
                                         Float.parseFloat(content[2]));
                    list.addZline(z);
                } else if (lastLine.charAt(0) == 'D') {
                    if(z != null) {
                        z.addDline(lastLine);
                    } else {
                        throw new RuntimeException("Wrong file format! Dline should always follow Zline");
                    }
               
                } else if (lastLine.charAt(0) == 'I') {
                    list.addIline(lastLine);
                } else if (lastLine.charAt(0) >= '0' && lastLine.charAt(0) <= '9') {
                    break; // digits, peaks starts
                } else {
                    continue;  // ignor lines other than I, D, S or peaks, e.g., empty line 
                }
            }
        } else {
            // process file in old format, i.e, spectra delimiter is ":"
           do {
                String line = lastLine.substring(1);
                content = line.split("\\.");

                //int lowscan = Integer.parseInt(content[1]);
                //int hiscan = Integer.parseInt(content[2]);

                int lowscan = Integer.parseInt(content[0]);
                int hiscan = Integer.parseInt(content[1]);
                //System.out.println("line: " + line + "\tlowscan: " + lowscan + "\thiscan: " + hiscan);

                list.setLoscan(lowscan);
                list.setHiscan(hiscan);
                if (fileFormat.equals("ms2")) {
  //                  list.setPrecursorMass(Float.parseFloat(content[3]));
                }
                 lastLine = br.readLine();
                 String[] oldZline = lastLine.split(" ");
                
                 int chargestate = Integer.parseInt(oldZline[1]);
                 float mplush = Float.parseFloat(oldZline[0]);

                 //System.out.println("line: " + line + "\tlowscan: " + lowscan + "\thiscan: " + hiscan + "\tcharge: " + chargestate + "\tmplush: " + mplush);
                 //z = new Zline(Integer.parseInt(oldZline[1]), Float.parseFloat(oldZline[0]));
                 z = new Zline(chargestate, mplush);
                 list.addZline(z);

                 //System.out.println("line: " + line + "\tlowscan: " + lowscan + "\thiscan: " + hiscan + "\tcharge: " + chargestate + "\tmplush: " + mplush +"\tnum zlines: " + list.getNumZlines());
            } while ((lastLine = br.readLine()) != null && lastLine.charAt(0) == spectraDelimiter);
        }
    }
}
