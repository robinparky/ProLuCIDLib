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
    private  List<String> accessionList = new ArrayList<>();
    private  List<String> descriptionList= new ArrayList<>();
    private TFloatArrayList mzList = new TFloatArrayList();
    private TFloatArrayList intensityList = new TFloatArrayList();
    private PeakList peakList = null;

    public LibrarySpectra(String sequence, int chargeState, float mz, float retTime, float score, float deltaCn,
                          String key, String filename, int scan, long id, String accession, String proteinDescription) {
        this.sequence = sequence;
        this.chargeState = chargeState;
        this.mz = mz;
        this.retTime = retTime;
        this.score = score;
        this.deltaCn = deltaCn;
        this.key = key;
        this.filename = filename;
        this.scan = scan;
        this.id = id;
        accessionList.add(accession);
        descriptionList.add(proteinDescription);
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
        accessionList.add(accession);
        descriptionList.add(description);
    }

    public List<String> getAccessionList() {
        return accessionList;
    }

    public List<String> getDescriptionList() {
        return descriptionList;
    }
}
