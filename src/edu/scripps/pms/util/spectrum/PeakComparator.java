/**
 * @file PeakComparater.java
 * This is the source file for edu.scripps.pms.util.spectrum.PeakComparater
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:40:54 $
 */



package edu.scripps.pms.util.spectrum;

import java.util.Comparator;
//import java.util.List;

public class PeakComparator implements Comparator {

    // how to comapare two PeakComparaters
    private static boolean compareByIntensity = false; 

    public PeakComparator(boolean compareByIntensity) {
        this.compareByIntensity = compareByIntensity;
    }   
    public void setCompareMode(boolean compareMode) {
        compareByIntensity = compareMode;
    }


    public int compare(Object peak1, Object peak2) {

        Peak p1 = (Peak) peak1;
        Peak p2 = (Peak) peak2;
        double f = 0;  // difference
        if (compareByIntensity) {
           f = p2.getIntensity() - p1.getIntensity();
        } else {
            f = p2.getM2z() - p1.getM2z();
        }
        
        if (f > 0) {
            return 1;
        } else if (f < 0) {
            return -1;
        } else {
            return 0;
        }
    }

    /**
     * Return true if both the intensities and m2zs are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {
        PeakComparator c = (PeakComparator)o;
        if(o == null) {
            return false;
        }
        return this.compareByIntensity == c.compareByIntensity;
    }
}
