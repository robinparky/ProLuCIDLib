/**
 * @file Peak.java
 * This is the source file for edu.scripps.pms.util.spectrum.Peak
 * @author Tao Xu
 * @date $Date: 2007/03/09 18:21:09 $
 */



package edu.scripps.pms.util.spectrum;

//import java.util.ArrayList;
//import java.util.List;

public class Peak {

    // how to comapare two Peaks

    private double m2z;
    private double intensity;
    private int chargeState = 0; // default value
    //private Header;
    // the order this Peak in the spectra, -1 indicates unknown 
    private int index = -1;
   
 
    public Peak(double m2z, double intensity) {
        this.m2z = m2z;
        this.intensity = intensity;
    }   
    public Peak(double m2z, double intensity, int chargeState) {
        this.m2z = m2z;
        this.intensity = intensity;
        this.chargeState = chargeState; // overwrite the default chargeState
    }

    public Peak(String line) {
        String [] arr = line.split(" ");
        m2z = Double.parseDouble(arr[0]);
        intensity = Double.parseDouble(arr[1]);
    }

    public int getChargeState() {
        return chargeState;
    }
    public double getM2z() {
        return m2z;
    }

    public double getIntensity() {
        return intensity;
    }
    // return the order this Peak in the PeakList. 
    // Return -1 is the order is unknown.
    // 0 means the first Peak in the list
    public int getIndex() {
        return index;
    }
    public void setM2z(double m2z) {
        this.m2z = m2z;
    }
    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }
    public void setChargeState(int chargeState) {
        this.chargeState = chargeState;
    }
    public void setIndex(int i) {
        index = i;
    }
    /**
     * Return true if both the intensities and m2zs are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {

        if (o == null) { 
            return false;
        }
        Peak p = (Peak)o;

        return this.m2z == p.m2z && this.intensity == p.intensity;
    }
}
