/**
 * @file MzxmlPeakList.java
 * This is the source file for edu.scripps.pms.util.spectrum.MzxmlPeakList
 * @author Tao Xu
 * @date $Date: 2011/01/29 00:12:53 $
 */



package edu.scripps.dia.util;

import edu.scripps.pms.util.spectrum.PeakList;

import java.util.ArrayList;
import java.util.Iterator;

public class MzxmlPeakList extends PeakList {
    private float totalIntensity = 0;
    private String listType;
    private String encodedM2zAndIntensities;
    private int encodingPrecision;
    private int index = 0;
    
    private ArrayList<MzxmlPeakList> peakLists = new ArrayList<MzxmlPeakList>();
    private int msLevel = 1; // for ms1
    private double retentionTime = 0;
    private double basePeakM2z = 0;
    private double basePeakIntensity = 0;
    private double totalIonCurrent = 0;

    private int precursorScanNum; 
    private double precursorIntensity; 

    public void addChildSpectrum(MzxmlPeakList mpl) {
        peakLists.add(mpl);
    }
    public int getNumChildSpectrum() {
        return peakLists.size();
    }
    public Iterator<MzxmlPeakList> getChildren() {
        return peakLists.iterator();
    }
    public void setEncodedM2zAndIntensities(String s, int precision) {
        encodedM2zAndIntensities = s;
        encodingPrecision = precision;
    }
    public String getEncodedM2zAndIntensities() {
        return  encodedM2zAndIntensities;
    }

    public int getEncodingPrecision() {
        return encodingPrecision;
    }
    public void setEncodingPrecision(int ep) {
        encodingPrecision = ep;
    }
    public double getBasePeakM2z() {
        return basePeakM2z;
    }
    public void setBasePeakM2z(double m2z) {

        basePeakM2z = m2z;
    }
    public void setMsLevel(int level) {
        msLevel = level;
    }  
    public int getMsLevel() {
        return msLevel;
    }  
    public double getBasePeakIntensity() {
        return basePeakIntensity;
    }
    public void setBasePeakIntensity(double d) {
        basePeakIntensity = d;
    }
    public void setTotalIonCurrent(double d) {
        totalIonCurrent = d;
    }
    public double getTotalIonCurrent() {
        return totalIonCurrent;
    }
    public void setRetentionTime(double d) {
        this.retentionTime = d;
    }
    public double getRetentionTime() {
        return retentionTime;
    }
    public void setPrecursorIntensity(double d) {
        precursorIntensity = d;
    } 
    public int getPrecursorScan() {
        return precursorScanNum; 
    }
    public void setPrecursorScan(int num) {
        precursorScanNum = num; 
    }


}
