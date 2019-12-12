/**
 * @file PeakList.java
 * This is the source file for edu.scripps.pms.util.spectrum.PeakList
 * @author Tao Xu
 * @date $Date: 2011/01/29 00:12:53 $
 */



package edu.scripps.pms.util.spectrum;

import edu.scripps.pms.mspid.MassSpecConstants;
import gnu.trove.TFloatArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.DoubleStream;


public class PeakList {
    public static final String RETENTIONTIME = "I\tRetTime";
    public static final String PRECURSORSCAN = "I\tPrecursorScan";
    public static final String PRECURSORINT = "I\tPrecursorInt"; // precursor intensity
    public static final String IONINJECTIONTIME = "I\tIonInjectionTime";
    public static final String INSTRUMENTTYPE = "I\tInstrumentType";
    public static int DEFAULTNUMPEAKS = 1000;
    public static final int DEFAULTSPECTRUMSIZE = 100000;
    public static boolean SORTBYM2Z = false;
    public static boolean SORTBYINTENSITY = true;
    private double totalIntensity = 0;
    private int loscan;
    private int hiscan;
    private double precursorMass; // this is actually precursor m/z, not M+H
    private List<String> hlines;
    private List<String> ilines = new ArrayList<String>();
    private ArrayList<Zline> zlines = new ArrayList<Zline>();
    private ArrayList<Peak> peaks = new ArrayList<Peak>(DEFAULTNUMPEAKS);
    private String listType;
    private int index = 0;
    private boolean isEtdSpectrum;
    private String activationType = null; // fragmentaion method CID or ETD
    private String instrumentType = null; // instrument type ITMS or FTMS
    private ArrayList<Peak> peaksSortedByIntensity = null;
    private ArrayList<Peak> peaksSortedByM2z = null;
    private double retentionTime = -1;
    private double ionInjectionTime = -1;
    private int precursorScanNumber = -1;
    private double precursorint = -1;
    private String filename;

    private double minM2z = Double.NaN;
    private double maxM2z = Double.NaN;

    public static final DecimalFormat threeDigits = new DecimalFormat("0.000");
    public static final DecimalFormat fourDigits = new DecimalFormat("0.0000");
    public static final DecimalFormat fiveDigits = new DecimalFormat("0.00000");
    public static final DecimalFormat twoDigits = new DecimalFormat("0.00");
    public static final DecimalFormat oneDigit = new DecimalFormat("0.0");

    // determine if the if this is a +1 spectrum. If so, keep only the +1 Zline, otherwise keep +2 and +3 Zlines
    public void reAssignCharge123() {
       if(getNumZlines() == 1) return;
       double allintensity = 0;
       double largemzintensity = 0;  // for intensity of ions with mz greater than precursor mz
       Iterator<Peak> it = peaks.iterator();
       while(it.hasNext()) {
           Peak p = it.next();
           if(p != null) {
               allintensity += p.getIntensity();
               if(p.getM2z() > (precursorMass + 1)) {
                   largemzintensity += p.getIntensity();
               }

           }
       }

       boolean ischarge1 = largemzintensity/allintensity < 0.05? true : false;

       ArrayList<Zline> temzlines = zlines;
       zlines = new ArrayList();
       Iterator<Zline> zit = temzlines.iterator();
       if(ischarge1) {
           while(zit.hasNext()) {
               Zline z = zit.next();
               if(z.getChargeState() == 1) {
                   zlines.add(z);
                   break;
               }
           }
       } else {
           while(zit.hasNext()) {
               Zline z = zit.next();
               int charge = z.getChargeState();
               if(charge == 2 || charge == 3 ) {
                   zlines.add(z);
               }
           }

       }

    }

    public int getNumZlines() {
        return zlines.size();
    }
    public Peak getPeak(double m2z, double massTolerance) {
        int index = getIndex(m2z, massTolerance);
        if(index > -1) {
            return getSortedPeaks(SORTBYM2Z).get(index);
        } else {
            return null;
        }
    }
    public boolean isEtdSpectrum() {
        return isEtdSpectrum;
    }
    public void resetZlines(ArrayList<Zline> zs) {
        zlines = zs;
    }

    // return -1 if the peak not fund
    private int getIndex(double m2z, double massTolerance) {
        ArrayList<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        int index = -1;
        int max = peaks.size() - 1;
        int min = 0;
        int mid = max/2;
        double maxM2z = m2z + massTolerance;
        double minM2z = m2z - massTolerance;

        while(max > (min+1)) {
            double m = sortedPeaks.get(mid).getM2z();
            if(m >= minM2z && m <= maxM2z) {
                return mid;
            }

            if(m > minM2z) {
                max = mid;
                mid = (min + mid)/2;
            } else if(m < minM2z) {
                min = mid;
                mid = (mid + max)/2;
            }
        }
        return index;
    }
    public double getIonInjectionTime() {
        if(ionInjectionTime == -1) {
            for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
                String i = it.next();
                if(i.startsWith(IONINJECTIONTIME)) {
                    String[] contents = i.split("\t");
                    ionInjectionTime = Double.parseDouble(contents[2]);
                    break;
                }

            }
        }
        return ionInjectionTime;
    }
    public String getInstrumentType() {

        if(instrumentType == null) {
            for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
                String i = it.next();
                if(i.startsWith(INSTRUMENTTYPE)) {
                    String[] contents = i.split("\t");
                    instrumentType = contents[2];
                    break;
                }

            }
        }
        return instrumentType;
    }
    public double getPrecursorInt() {
        if(precursorint == -1) {
            for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
                String i = it.next();
                if(i.startsWith(PRECURSORINT)) {
                    String[] contents = i.split("\t");
                    precursorint = Double.parseDouble(contents[2]);
                    break;
                }

            }
        }
        return precursorint;
    }
    public int getPrecursorScanNumber() {
        if(precursorScanNumber == -1) {
            for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
                String i = it.next();
                if(i.startsWith(PRECURSORSCAN)) {
                    String[] contents = i.split("\t");
                    precursorScanNumber = Integer.parseInt(contents[2]);
                    break;
                }

            }
        }
        return precursorScanNumber;
    }
    public double getRetentionTime() {
        if(retentionTime == -1) {
            for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
                String i = it.next();
                if(i.startsWith(RETENTIONTIME)) {
                    String[] contents = i.split("\t");
                    retentionTime = Double.parseDouble(contents[2]);
                    break;
                }
            }
        }
        return retentionTime;
    }
    public void setHlines(List<String> hlines) {
       this.hlines = hlines;
    }

    public Peak getPeakByM2z(double m2z) {
        double minDiff = 1000;
        Peak prcPeak = null;
        for(Iterator<Peak> it = peaks.iterator(); it.hasNext();) {
            Peak currPeak = it.next();
            double m2zCurrent = currPeak.getM2z();
            double diff = m2z > m2zCurrent? m2z-m2zCurrent : m2zCurrent-m2z; // avoid abs()
            if(diff == 0) {
                return currPeak;
            }
            if(diff < minDiff) {
                minDiff = diff;
                prcPeak = currPeak;
            }
        }
        return prcPeak;
    }
    public int getPrecursorScan() {
        for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
            String s = it.next();
            if(s.startsWith("I\tPrecursorScan")) {
                String[] elements = s.split("\t");
                return Integer.parseInt(elements[2]);
            }
        }
        return -1; // I line for PrecursorScan was not found
    }
    public void setZlines(ArrayList<Zline> zlines) {
       this.zlines = zlines;
    }
    public void setLoscan(int loscan) {
        this.loscan = loscan;
    }
    public void setHiscan(int hiscan) {
        this.hiscan = hiscan;
    }
    public List<String> getIlines() {
        return ilines;
    }
    public List<String> getHlines() {
        return hlines;
    }
    public void setPrecursorMass(double precursorMass) {
        this.precursorMass = precursorMass;
    }
    public String getActivationType() {
        return activationType;
    }
    public void addIline(String l) {
        if(l == null) {
            return;
        }
        String arr [] = l.split("\t");
        //if(arr.length > 2 && arr[2].startsWith("ETD")) {
        if(arr.length > 2 && arr[1].startsWith("ActivationType")) {
            activationType = arr[2];
            if(activationType.startsWith("ETD") || activationType.startsWith("etd")) {
                isEtdSpectrum = true;
            }
//System.out.println("isETD: " + isEtdSpectrum);
        } else {
  //System.out.println("in addIline: " + l + "\tisETD: " + isEtdSpectrum);
        }
        ilines.add(l);
    }
    public void addZline(Zline z) {
        zlines.add(z);
    }

    public void addPeak(Peak p) {
        p.setIndex(index++);
        peaks.add(p);
        totalIntensity += p.getIntensity();
    }

    public int numPeaks() {
        return peaks.size();
    }
    public double getTotalIntensity() {
        return totalIntensity;
    }
    public Iterator<Zline> getZlines() {
        return zlines.listIterator();
    }
    public void sortPeaks(boolean sortByIntensity) {
        peaks = getSortedPeaks(sortByIntensity);
    }

    // group splitted peaks and retain 10 peak for every 100 m/z, assume the spectrum is sorted by m/z
    public void processHcdSpectrum() {
//System.out.println("Number of peaks before: " + peaks.size());
        if(peaks.size () < 10) return;

        ArrayList<Peak> newpeaks = new ArrayList<Peak>(DEFAULTNUMPEAKS);
        ArrayList<Peak> toberemoved = new ArrayList();

        Peak prev = peaks.get(0);
        for(int i = 1; i < peaks.size(); i++) {
            Peak p = peaks.get(i);
            if(p.getM2z() - prev.getM2z() < 0.04) {
                if(p.getIntensity() < prev.getIntensity()) {
                    toberemoved.add(p);
                } else {
                    toberemoved.add(prev);
                }
            }
            prev = p;
        }

        Iterator<Peak> it = toberemoved.iterator();
        while(it.hasNext()) {
            peaks.remove(it.next());
        }

        double maxmz = getMaxM2z();
        int maxpeakindex = peaks.size();
        int peakindex = 0;
        for(int i = 0; i*100 < maxmz; i++) {
            ArrayList<Peak> temp = new ArrayList(100);
            double uppermzlimit = (i+1)*100;
            while(peakindex < maxpeakindex) {
                Peak p = peaks.get(peakindex++);

                if(p.getM2z() <  uppermzlimit) {

                    temp.add(p);
                } else {
                    break;
                }
            }

//System.out.println("New Peak List:");
            //sort temp by intensity
            Collections.sort(temp, new PeakComparator(true));
            for(int j = temp.size() - 18; j >=0 && j < temp.size(); j++) {
                newpeaks.add(temp.get(j));
//System.out.println(temp.get(j).getM2z() + " " + temp.get(j).getIntensity());
            }

        }
//System.out.println("Number of peaks before: " + peaks.size());
        peaks = newpeaks;
//System.out.println("Number of peaks after: " + peaks.size());
        peaksSortedByIntensity = null;
        peaksSortedByM2z = null;
    }

    /**
     * Return a sorted list (incremental by m2z or intensity)
     * for the Peaks in this PeakList
     * @param sortByIntensity - indicate how the list should be sorted,
     *                          true for sort by intensity,
     *                          false for sort by M2z
     * @note user must not modify the List returned.
     *
     */
    public synchronized ArrayList<Peak> getSortedPeaks(boolean sortByIntensity) {
        ArrayList<Peak> sortedPeaks =  new ArrayList(DEFAULTNUMPEAKS);
        if (sortByIntensity) {
            if (peaksSortedByIntensity != null) {
                return peaksSortedByIntensity;
            } else {
                peaksSortedByIntensity = sortedPeaks;
            }
        } else {
            if (peaksSortedByM2z != null) {
                return peaksSortedByM2z;
            } else {
                peaksSortedByM2z = sortedPeaks;
            }
        }

        for (Peak p : peaks) {
            sortedPeaks.add(p);
        }

        Collections.sort(sortedPeaks, new PeakComparator(sortByIntensity));
        return sortedPeaks;
    }
    public synchronized ArrayList<Peak> getSortedPeaks(int numPeaks, boolean sortByIntensity) {
        ArrayList topPeaks = new ArrayList<Peak>(numPeaks);

        List<Peak> sortedList = getSortedPeaks(sortByIntensity);
        int totalPeaks = sortedList.size()-1;
        while(numPeaks > 0 && totalPeaks >= 0) {
            Peak p = sortedList.get(totalPeaks--);
            //if(p.getM2z() > getMaxChargeState()*getPrecursorMass()){
            topPeaks.add(p);
            numPeaks--;

        }
        return topPeaks;
    }
    public ListIterator<Peak> getPeaks() {
        return peaks.listIterator();
    }

    public String getListType()
    {
        return listType;
    }

    public void setListType(String listType)
    {
        this.listType = listType;
    }
    /**
     * Return the low scan number of this
     */
    public int getLoscan() {
        return loscan;
    }

    public int getHiscan() {
        return hiscan;
    }

    // return the m/z of the precursor ion, not the M+H value
    public double getPrecursorMass() {
        return precursorMass;
    }
    public double getPrecursorMz() {
        return precursorMass;
    }
    public double getMaxPrecursorMass() {
        double maxPrecMass = 0;
        for(Zline z : zlines) {
            double precMass = z.getM2z();
            maxPrecMass = maxPrecMass > precMass? maxPrecMass : precMass;
        }
        return maxPrecMass;
    }
    public double getMaxM2z() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        return sortedPeaks.get(0).getM2z();
    }
    public double getMinM2z() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        return sortedPeaks.get(numPeaks()-1).getM2z();
    }
    public double getMaxIntensity() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
        return sortedPeaks.get(sortedPeaks.size()-1).getIntensity();
    }
    public double getMinIntensity() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
        return sortedPeaks.get(0).getIntensity();
    }
    // get the intensity mean of the percentage percent of least intense peak, usuall can be used as noise level
    public double getNoiseLevel(int percentage) {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        int numPeaks = (int)(peaks.size()*percentage/100 + 0.5);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks + "\tpercentage: " + percentage);
        numPeaks = numPeaks > 0? numPeaks : 1;
        double sumnoiseint = 0;
        int lastIndex = peaks.size() - 1;
        for(int i = 0; i < numPeaks; i++) {
//System.out.println(sortedPeaks.get(i).getIntensity());
            //leastIntensePeaks.add(sortedPeaks.get(i));
            sumnoiseint += sortedPeaks.get(i).getIntensity();
        }
        return sumnoiseint/numPeaks;
    }
    public ArrayList<Peak> getLeastIntensePeaks(int percentage) {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        int numPeaks = (int)(peaks.size()*percentage/100 + 0.5);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        numPeaks = numPeaks > 0? numPeaks : 1;
        ArrayList<Peak> leastIntensePeaks = new ArrayList<Peak>(numPeaks);
        int lastIndex = peaks.size() - 1;
        for(int i = 0; i < numPeaks; i++) {
//System.out.println(sortedPeaks.get(i).getIntensity());
            leastIntensePeaks.add(sortedPeaks.get(i));
        }
        return leastIntensePeaks;
    }
    // assume the peaks are sorted by m2z
    public double getMaxIntensity(double minM2z, double maxM2z) {
        double maxIntensity = -10000;
        for(Iterator<Peak> it = peaks.iterator(); it.hasNext();) {
            Peak p = it.next();
            double intensity = p.getIntensity();
            double m2z = p.getM2z();
            if(m2z >= minM2z) {
                if(m2z <= maxM2z) {
                    maxIntensity = maxIntensity > intensity? maxIntensity : intensity;
                } else {
                    break;
                }
            }
        }
        return maxIntensity;
    }

    public Peak getMinimumDeltaMassPeak(double mz, double tol) {
        int index = getIndex(mz, tol);
        if(index == -1) {
            return null;
        }
        double maxM2z = mz + tol;
        double minM2z = mz - tol;
        ArrayList<Peak> ps = getSortedPeaks(SORTBYM2Z);
        int maxindex = ps.size();
        Peak p = ps.get(index);
        Peak result = p;
        Peak tempp = p;
        int tempindex = index;
        double bestdeltamass = Math.abs(p.getM2z() - mz);
        while(tempp.getM2z() < maxM2z && tempindex < maxindex) {
            tempp = ps.get(++tempindex);
            double deltamass = Math.abs(tempp.getM2z() - mz);
            if(deltamass <= bestdeltamass) {
                bestdeltamass = deltamass;
                result = tempp;
            } else {
                break;
            }
        }
        tempp = p;
        tempindex = index;
        while(tempp.getM2z() > minM2z && tempindex > 0) {
            tempp = ps.get(--tempindex);
            double deltamass = Math.abs(tempp.getM2z() - mz);
            if(deltamass <= bestdeltamass) {
                bestdeltamass = deltamass;
                result = tempp;
            } else {
                break;
            }

        }
        return result;
    }
    // assume the peaks are sorted by m2z, use charge state to store numpeaksinrange
    public Peak getMaxIntensePeak(double minM2z, double maxM2z) {
        double maxIntensity = -10000;
        Peak maxIntensePeak = null;
        int numpeaksinrange = 0;
        for(Iterator<Peak> it = peaks.iterator(); it.hasNext();) {
            Peak p = it.next();
            double intensity = p.getIntensity();
            double m2z = p.getM2z();
            if(m2z >= minM2z) {
                if(m2z <= maxM2z) {
                    numpeaksinrange++;
                    if(maxIntensity < intensity) {
                        maxIntensity = intensity;
                        maxIntensePeak = p;
                    }
                } else {
                    break;
                }
            }
        }
        if(maxIntensePeak != null) {
            maxIntensePeak.setChargeState(numpeaksinrange);
        }
        return maxIntensePeak;
    }
    // mainly for removing precursor ions of ETD data
    public ArrayList<Range> getPrecursorRanges(Zline z) {
        ArrayList<Range> precs = new ArrayList<Range>();
        int chargestate = z.getChargeState();
        double mass = (precursorMass - MassSpecConstants.MASSH)*chargestate;
        if(isEtdSpectrum) {
            for(int i = chargestate; i > 0; i--) {
                double mz = (mass +  MassSpecConstants.MASSH*i)/i;
                double low = mz - 35/i - 2;
                precs.add(new Range(low, mz+2));
            }
        } else {
            precs.add(new Range(precursorMass-(35/chargestate)-2, precursorMass+2));
        }
        return precs;
    }
    public double getAvgIntensityOfLeastIntensePeaks(int percentage) {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        int numPeaks = (int)(peaks.size()*percentage/100 + 0.5);

        numPeaks = numPeaks > 0? numPeaks : 1;
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        //int lastIndex = peaks.size() - 1;
        double totalIntens = 0;
        for(int i = 0; i < numPeaks; i++) {
            totalIntens += sortedPeaks.get(i).getIntensity();
//System.out.println(i + "\t" + sortedPeaks.get(i).getIntensity() +"\ttotalIntens: " + totalIntens + "\tnumPeaks: " + numPeaks);
        }
//System.out.println("totalIntensity: " + totalIntens + "\tumPeaks: " + numPeaks + "\tavgNoise: " + totalIntens/numPeaks);
        return totalIntens/numPeaks;
    }
    public String getSpectrumWithHlines() {
        StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
        for(String h : hlines) {
            sb.append(h);
            sb.append("\n");
        }
        getSpectrumWithoutHlines(sb);
        return sb.toString();
    }
    public String getSpectrumWithoutHlines() {
        StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
        getSpectrumWithoutHlines(sb);
        return sb.toString();
    }
    // everything except the peaks
    public void getScanInfo(StringBuffer sb) {
        addSline(sb);
        addIlines(sb);
        addZlines(sb);
    }
    private void getSpectrumWithoutHlines(StringBuffer sb) {
        addSline(sb);
        addIlines(sb);
        addZlines(sb);
        addSpectrum(sb);
    }
    private void addIlines(StringBuffer sb) {
        //for(String i : ilines) {
        for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
            String i = it.next();
            sb.append(i);
            sb.append("\n");
        }
    }
    private void addSline(StringBuffer sb) {
//System.out.println("loscan: " + loscan);
        sb.append("S");
        sb.append("\t");
        sb.append(loscan);
        sb.append("\t");
        sb.append(hiscan);
        if(getNumZlines() != 0) {
            sb.append("\t");
            sb.append(fiveDigits.format(precursorMass));
        }
        sb.append("\n");
    }
    private void addZlines(StringBuffer sb) {

        //for(Zline z : zlines) {
        for(Iterator<Zline> it = zlines.iterator(); it.hasNext();) {
            Zline z = it.next();
            sb.append("Z");
            sb.append("\t");
            sb.append(z.getChargeState());
            sb.append("\t");
            sb.append(fiveDigits.format(z.getM2z()));
            sb.append("\n");
            for(String d : z.getDlines()) {
                sb.append(d);
                sb.append("\n");
            }
        }
    }
    private void addSpectrum(StringBuffer sb) {
        //for(Peak p : peaks) {
        for(Iterator<Peak> it = peaks.iterator(); it.hasNext();) {
            Peak p = it.next();
            sb.append(fiveDigits.format(p.getM2z()));
            sb.append(" ");
            sb.append(oneDigit.format(p.getIntensity()));
            sb.append("\n");
        }
//System.out.println(sb.toString());
    }
    public String getFragmentIons() {
        StringBuffer sb = new StringBuffer(2000);
        //for(Peak p : peaks) {
        for(Iterator<Peak> it = peaks.iterator(); it.hasNext();) {
            Peak p = it.next();
            sb.append(fiveDigits.format(p.getM2z()));
            sb.append(" ");
            sb.append(oneDigit.format(p.getIntensity()));
            sb.append("\n");
        }
        return sb.toString();
    }
    public int getMaxChargeState() {
        int maxChargeState = 0;
        for(Zline z : zlines) {
            if(z.getChargeState() > maxChargeState) {
                maxChargeState = z.getChargeState();
            }
        }
        return maxChargeState;
    }
    public int getFirstChargeState() {
        for(Zline z : zlines) {
            return z.getChargeState();
        }
        return 2; // default 
    }
    private List<Peak> getMostIntensPeaks(int numPeaks, double lowM2zLimit, double highM2zLimit) {
        ArrayList<Peak> topPeaks = new ArrayList<Peak>(numPeaks);
        List<Peak> sortedList = getSortedPeaks(true);
        int totalPeaks = sortedList.size()-1;
        while(numPeaks > 0 && totalPeaks >= 0) {
            Peak p = sortedList.get(totalPeaks--);
            //if(p.getM2z() > getMaxChargeState()*getPrecursorMass()){
            double m2z = p.getM2z();
            if(m2z > lowM2zLimit && m2z < highM2zLimit){
                topPeaks.add(p);
                numPeaks--;
            }

        }
        return topPeaks;
    }

    public PointList [] calcQCorrs(int numPeaks, int accuracyFactor, double massWindow) {
        int massH = (int)(MassSpecConstants.MASSH*accuracyFactor + 0.5);
        PointList [] points = new PointList[getNumZlines()];
        int zlineCounter = 0;

        //for(Zline z : zlines) {
        for(Iterator<Zline> it = zlines.iterator(); it.hasNext();) {
            Zline z = it.next();
            int chargeState = z.getChargeState();
            double lowLimit = massWindow*chargeState;
            double noShiftMass = precursorMass*chargeState;
            double highLimit = noShiftMass - lowLimit;
            List<Peak> topList = getMostIntensPeaks(numPeaks, lowLimit, highLimit);
            int numBins = (int)(precursorMass*chargeState*accuracyFactor + 0.5);
            double [] intens = new double[numBins];
            for(Peak p : topList) {
                int index = (int)(p.getM2z()*accuracyFactor+0.5);
                intens[index] += p.getIntensity();
            }
            //double [] temp = intens;
            intens = reverse(intens); // get the revered list for correlation
            int numShift = (int)massWindow*chargeState*accuracyFactor;
            int offSet = numShift/2;
//System.out.println("scan#: " + getLoscan() +"\tprecuMass: " + precursorMass +  "\toffSet: " + offSet);
            double [] qcorr = new double [numShift];

            for(Peak p : topList) {
                int m2z = (int)(p.getM2z()*accuracyFactor + 0.5);
                double intensity = p.getIntensity();
                int startIndex = m2z-offSet;
                for(int i = 0; i < numShift; i++) {
                    qcorr[i] += intensity*intens[startIndex + i];
                }
            }
            PointList pointList = new PointList();
            double firstMass = noShiftMass+massWindow/2*chargeState;
            double dAccFactor = accuracyFactor; // convert to accuracyFactor to double
System.out.println("FirstMass: " + firstMass + "\tnoShiftMass: " + noShiftMass+ "lowLimit: " + lowLimit + "\tMax: " + (numBins+offSet)/dAccFactor + "\tMin: " + (numBins-offSet)/dAccFactor);
            for(int i = 0; i < qcorr.length; i++) {
                pointList.addPoint(new Point(firstMass-i/dAccFactor, qcorr[i]));
            }
            points[zlineCounter++] = pointList;
            int indexMax = getMaxIndex(qcorr);
            //double mass = (numBins-offSet-indexMax-1.5*massH)/(accuracyFactor+0.0); 
            double mass = firstMass-indexMax/(accuracyFactor+0.0);
//System.out.println("chargeState: " + z.getChargeState()+"\tnumBins: " + numBins +"\tm2z: " + z.getM2z() + "\tmass: " + mass + "\tindexMax: " + indexMax + "\tqcorr: " + qcorr[indexMax]);
//System.out.println(z.getChargeState() + "\t"+qcorr[indexMax] + "\t" + (z.getM2z()-mass));
            //String dline = "D\t" + mass + "\t" + qcorr[indexMax];
            //z.addDline(dline);
        }
        return points;
    }
    private double [] reverse(double [] intens) {
        double [] reversed = new double[intens.length];
        int lastIndex = intens.length - 1;
        for(int i = 0; i <= lastIndex; i++) {
            reversed[lastIndex-i] = intens[i];
        }
        return reversed;
    }
    private int getMaxIndex(double [] qcorr) {

        int indexMax = 0;
        double max = 0;
        for(int i = 0; i < qcorr.length; i++) {
            if(qcorr[i] > max) {
                max = qcorr[i];
                indexMax = i;
            }

        }
//System.out.println("max: " + max);
        return indexMax;
    }

    public void setActivationType(String activationType) {

        this.activationType = activationType;
        if(activationType.startsWith("ETD") || activationType.startsWith("etd")) {
            isEtdSpectrum = true;
        }
    }

    public void setFilename(String filename)
    {
        this.filename = filename;
    }

    public String getFilename()
    {
        return filename;
    }

    public TFloatArrayList generateBinSpectra(int shift, float fragment_tolerance,
                                              int intensities_sum_or_mean_or_median, float minMz, float maxMz) {
        TFloatArrayList bin_spec_al = new TFloatArrayList();

        float binSize = (fragment_tolerance * 2),
                upperLimit = maxMz + 0.00001f;
        for (float lowerLimit = minMz; lowerLimit < upperLimit; lowerLimit = lowerLimit + binSize) {
            float tmp_intensity_bin = 0;
            DescriptiveStatistics obj = new DescriptiveStatistics();
            for (Iterator<Peak> peakIterator =  this.getPeaks(); peakIterator.hasNext();) {
                Peak p = peakIterator.next();
                double mz = p.getM2z()+ shift;
                if (mz >= lowerLimit && mz < lowerLimit + binSize) {
                    obj.addValue(p.getIntensity());
                }
            }
            if (obj.getN() > 0) {
                if (intensities_sum_or_mean_or_median == 0) {
                    tmp_intensity_bin = (float)obj.getSum();
                } else if (intensities_sum_or_mean_or_median == 1) {
                    tmp_intensity_bin =(float) obj.getMean();
                } else if (intensities_sum_or_mean_or_median == 2) {
                    tmp_intensity_bin = (float)obj.getPercentile(50);
                }
            }
            // put every bin_pectrum
            bin_spec_al.add(tmp_intensity_bin);
        }
        // convert an arraylist to double array
        // initiate size of array

        return bin_spec_al;
    }

    public List<Peak> getPeaksList()
    {
        return peaks;
    }
}
