/**
 * @file PointComparater.java
 * This is the source file for edu.scripps.pms.util.spectrum.PointComparater
 * @author Tao Xu
 * @date $Date: 2006/01/27 23:46:05 $
 */



package edu.scripps.pms.util.spectrum;

import java.util.Comparator;
//import java.util.List;

public class PointComparator implements Comparator {

    // how to comapare two PointComparaters
    private static boolean compareByIntensity = false; 

    public PointComparator(boolean compareByIntensity) {
        this.compareByIntensity = compareByIntensity;
    }   
    public void setCompareMode(boolean compareMode) {
        compareByIntensity = compareMode;
    }


    public int compare(Object point1, Object point2) {

        Point p1 = (Point) point1;
        Point p2 = (Point) point2;
        double f = 0;  // difference
         
        if (compareByIntensity) {
           f = p1.getIntensity() - p2.getIntensity();
        } else {
            f = p1.getXValue() - p2.getXValue();
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
     * Return true if both the intensities and xValues are the same, 
     * otherwise return false
     */
    public boolean equals(Object o) {
        PointComparator c = (PointComparator)o;
        if(o == null) {
            return false;
        }
        return this.compareByIntensity == c.compareByIntensity;
    }
}
