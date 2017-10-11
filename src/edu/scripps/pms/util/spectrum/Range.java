

/**
 * @file Range.java
 * This is the source file for edu.scripps.pms.util.spectrum.Range
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.util.spectrum;

public class Range {
    private double lowBound;
    private double highBound;
 
    public Range(double low, double high) {
        lowBound = low;
        highBound = high;
    }
    public void setLowBound(double low) {
        lowBound = low;
    } 
    public void setHighBound(double high) {
        highBound = high;
    }
    public double getLowBound() {
        return lowBound;
    } 
    public double getHighBound() {
        return highBound;
    }
    public boolean isInRange(double value) {
        return value >= lowBound && value <= highBound;
    }
}



