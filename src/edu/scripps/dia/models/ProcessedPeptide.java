package edu.scripps.dia.models;

import gnu.trove.TIntArrayList;
import gnu.trove.TIntHashSet;

import java.util.*;

/**
 * Created by yateslab on 7/19/17.
 */
public class ProcessedPeptide {
    private String DatabasePath;
    private String sequence;
    private int chargeState;
    private String key;
    private Set<String> fileScanSet = new HashSet<>();
    private Map<String,TIntHashSet> fileScanMap = new HashMap<>();
    private Set<String> fileSet = new HashSet<>();
    private Set<String> proteinInfo = new HashSet<>();
    private float theoreticalMass;
    public ProcessedPeptide(String sequence, int chargeState)
    {
        init(sequence,chargeState);
    }

    public ProcessedPeptide(String sequence, float theoreticalmass, int chargeState, String file, int scan)
    {
        addFileScan(file,scan);
        this.theoreticalMass = theoreticalmass;
        init(sequence,chargeState);
    }

    private void init(String sequence, int chargeState)
    {
        this.sequence = sequence;
        this.chargeState = chargeState;
        key = sequence +chargeState;
    }

    public void addFileScan(String file, int scan)
    {
        TIntHashSet scanList = fileScanMap.get(file);
        if(scanList ==null)
        {
            scanList = new TIntHashSet();
            scanList.add(scan);
            fileScanMap.put(file,scanList);
        }
        else
        {
            scanList.add(scan);
        }
    }

    public void addProtienInfo(String locus, String description)
    {
        proteinInfo.add(locus+"\t"+description);
    }

    public Iterator<String> proteinIterator()
    {
        return proteinInfo.iterator();
    }



    @Override
    public int hashCode() {
        return key.hashCode();
    }


    public String getDatabasePath() {
        return DatabasePath;
    }

    public String getSequence() {
        return sequence;
    }

    public int getChargeState() {
        return chargeState;
    }

    public String getKey() {
        return key;
    }

    public Set<String> getFileScanSet() {
        return fileScanSet;
    }

    public Set<String> getProteinInfo() {
        return proteinInfo;
    }

    public Set<String> getFileSet() {
        return fileSet;
    }

    public Map<String, TIntHashSet> getFileScanMap() {
        return fileScanMap;
    }

    public float getTheoreticalMass() {
        return theoreticalMass;
    }
}
