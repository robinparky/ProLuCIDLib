package edu.scripps.dia;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import gnu.trove.TFloatArrayList;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by Titus Jung titusj@scripps.edu on 3/29/18.
 */
public class LibrarySpectra {

    public final String sequence;
    public final int chargeState;
    public final float mz;
    public final float retTime;
    public final float score;
    public final float deltaCn;
    public final String key;
    public final String filename;
    public final int scan;
    public final long id;
    public final int massKey;
    public final boolean isHeavy;
    private float startTime;
    private float endTime;

    private  List<String> accessionList = new ArrayList<>();
    private  List<String> descriptionList= new ArrayList<>();
    private TFloatArrayList mzList = new TFloatArrayList();
    private TFloatArrayList intensityList = new TFloatArrayList();
    private  PeakList peakList = null;
    private boolean isDecoy = false;


    public static String ReveseSequence(String sequence)
    {
        StringBuilder sb = new StringBuilder();
        sb.append(sequence.charAt(0)).append(".");
        int lastCutLoc = sequence.lastIndexOf(".")-1;
        while(!Character.isAlphabetic(sequence.charAt(lastCutLoc)))
        {
            lastCutLoc--;
        }
        String end = sequence.substring(lastCutLoc);
        String mid = sequence.substring(sequence.indexOf(".")+1, lastCutLoc );
        StringBuilder midSb = new StringBuilder();
        int loc = mid.indexOf("(");
        String modStr="";
        if(loc>-1)
            modStr = mid.substring(mid.indexOf("("), mid.lastIndexOf(")")+1);

        String cleanSeq = mid.replaceAll("\\([0-9.\\-]*\\)", "");
        if(cleanSeq.contains("("))
        {
            System.out.println();
        }
        for(int i=cleanSeq.length()-1; i>=0; i--)
        {
            midSb.append(cleanSeq.charAt(i));
        }
        if(loc>0)
        midSb.insert(loc,modStr);

        sb.append(midSb);
       // sb.append(".");
        sb.append(end);
        return  sb.toString();
    }




    public LibrarySpectra(String sequence, int chargeState, float mz, float retTime, float score, float deltaCn,
                          String key, String filename, int scan, long id, String accession, String proteinDescription,
                          int massKey) {
        this.sequence = sequence;
        this.chargeState = chargeState;
        this.mz = mz;
        this.retTime = retTime;
        this.score = score;
        this.deltaCn = deltaCn;
        this.key = key;
        this.isHeavy = key.charAt(key.length()-1) == 'H' ;
        this.filename = filename;
        this.scan = scan;
        this.id = id;
        this.massKey = massKey;
        accessionList.add(accession);
        descriptionList.add(proteinDescription);
       // this.accession = accession;
       // this.proteinDescription = proteinDescription;
    }

    public LibrarySpectra(int massKey, int cs,String sequence, boolean isDecoy, String key, float retTime ) {
        if(isDecoy)
        {
            sequence = ReveseSequence(sequence);
        }
        this.sequence = sequence;
        this.retTime =retTime;
        this.chargeState = cs;
        this.mz = massKey/1000.0f;

        this.score = 0;
        this.deltaCn = 0;
        this.key = key;
        this.isHeavy = key.charAt(key.length()-1) == 'H' ;
        this.filename = "";
        this.scan = 0;
        this.id = 0;
        this.massKey = massKey;
        this.isDecoy = isDecoy;

        // this.accession = accession;
        // this.proteinDescription = proteinDescription;
    }




    public TFloatArrayList getMzList() {
        return mzList;
    }

    public void setMzList(TFloatArrayList mzList) {
        this.mzList = mzList;
    }

    public TFloatArrayList getIntensityList() {
        return intensityList;
    }

    public void setIntensityList(TFloatArrayList intensityList) {
        this.intensityList = intensityList;
    }
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(sequence).append("\t")
                .append(chargeState).append("\t")
                .append(filename).append("\t")
                .append(scan).append("\t")
                .append(mz).append("\t")
                .append(retTime).append("\t")
                .append(score).append("\t")
                .append(deltaCn).append("\t");

        return sb.toString();
    }

    public String printSpectra()
    {
        StringBuilder sb = new StringBuilder();
        for(int i=0; i<mzList.size(); i++)
        {
            sb.append(mzList.get(i)).append("\t").append(intensityList.get(i)).append("\n");
        }
        return sb.toString();
    }

    public void createPeakList()
    {
        peakList = new PeakList();
        for(int i=0; i<mzList.size(); i++)
        {
            peakList.addPeak(new Peak(mzList.get(i),intensityList.get(i)));
        }
    }
    public PeakList createPeakListCopy()
    {
        PeakList peakList = new PeakList();
        for(int i=0; i<mzList.size(); i++)
        {
            peakList.addPeak(new Peak(mzList.get(i),intensityList.get(i)));
        }
        return peakList;
    }

    public synchronized List<Peak> getIntensPeaks(int num) {
        if(peakList==null) createPeakList();
        List<Peak> sortedPeaks = peakList.getSortedPeaks(PeakList.SORTBYINTENSITY);
        List<Peak> intensPeaks = new LinkedList<Peak>();
        int numPeaks = sortedPeaks.size();
        int numIntensPeaks = 0;
        while(numIntensPeaks < num && numIntensPeaks < numPeaks) {
            numIntensPeaks++;
            Peak p = sortedPeaks.get(numPeaks-numIntensPeaks);
            intensPeaks.add(p);
        }
        return intensPeaks;
    }

    public void addProtein(String accession, String description)
    {
        if(isDecoy)
        {
            accession= "Reverse_"+accession;
        }
        accessionList.add(accession);
        descriptionList.add(description);
    }

    public List<String> getAccessionList() {
        return accessionList;
    }

    public List<String> getDescriptionList() {
        return descriptionList;
    }


    public PeakList getPeakList() {
        if(peakList==null) createPeakList();

        return peakList;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public void setDecoy(boolean decoy) {
        isDecoy = decoy;
    }


    public float getStartTime() {
        return startTime;
    }

    public void setStartTime(float startTime) {
        this.startTime = startTime;
    }

    public float getEndTime() {
        return endTime;
    }

    public void setEndTime(float endTime) {
        this.endTime = endTime;
    }
}
