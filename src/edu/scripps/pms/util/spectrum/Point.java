/**
 * @file Point.java
 * This is the source file for edu.scripps.pms.util.spectrum.Point
 * @author John Venable 
 * @date $Date: 2006/01/27 23:46:05 $
 */



package edu.scripps.pms.util.spectrum;

//import java.util.ArrayList;
//import java.util.List;

public class Point {

    // how to comapare two Peaks

    private double xValue;
    private double intensity;
    
    //private Header;
    // the order this Peak in the spectra, -1 indicates unknown 
    private int index = -1;
   
    public Point(double xValue, double intensity) {
        this.xValue = xValue;
        this.intensity = intensity;
    }   
    public double getXValue() {
        return xValue;
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
    public void setXValue(double xValue) {
        this.xValue = xValue;
    }
    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }
    public void setIndex(int i) {
        index = i;
    }
    /**
     * Return true if both the intensities and xValuess are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {

        if (o == null) { 
            return false;
        }
        Point p = (Point)o;

        return this.xValue == p.xValue && this.intensity == p.intensity;
    }
}
