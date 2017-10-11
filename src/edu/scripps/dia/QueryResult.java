package edu.scripps.dia;

import gnu.trove.TFloatArrayList;
import gnu.trove.TIntArrayList;

/**
 * Created by yateslab on 8/17/17.
 */
public class QueryResult {

    public final String sequence;
    public final float mass;
    private TFloatArrayList massList;
    private TIntArrayList intensityList;
    private int [] massArray;

    public QueryResult(String seq, float mass)
    {
        this.sequence = seq;
        this.mass = mass;
    }

    public QueryResult(String seq, float mass, TFloatArrayList massList, TIntArrayList intensityList)
    {
        this.sequence = seq;
        this.mass = mass;
        this.intensityList = intensityList;
        this.massList = massList;
    }

    public QueryResult(String seq, float mass, int[] massArray)
    {
        this.sequence = seq;
        this.mass = mass;
        this.massArray = massArray;
    }

    public TFloatArrayList getMassList() {
        return massList;
    }

    public void setMassList(TFloatArrayList massList) {
        this.massList = massList;
    }

    public TIntArrayList getIntensityList() {
        return intensityList;
    }

    public void setIntensityList(TIntArrayList intensityList) {
        this.intensityList = intensityList;
    }

    public int[] getMassArray() {
        return massArray;
    }

    public void setMassArray(int[] massArray) {
        this.massArray = massArray;
    }
}
