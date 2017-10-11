/**
 * @file PointList.java
 * This is the source file for edu.scripps.pms.util.spectrum.PointList
 * @author Tao Xu
 * @date $Date: 2006/01/27 23:46:05 $
 */



package edu.scripps.pms.util.spectrum;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.ListIterator;

public class PointList {

    public static int DEFAULTNUMPEAKS = 1000;
    public static boolean SORTBYXVALUE = false;
    public static boolean SORTBYINTENSITY = true;
   
    private double totalIntensity = 0;
    private ArrayList<Point> points = new ArrayList(DEFAULTNUMPEAKS);
    private String listType;
    private int index = 0;

    private ArrayList<Point> pointsSortedByIntensity = null;
    private ArrayList<Point> pointsSortedByXValue = null;

    public void addPoint(Point p) {
        p.setIndex(index++);
        points.add(p);
        totalIntensity += p.getIntensity();
    }

    public int numPoints() {
        return points.size();
    }
    public double getTotalIntensity() {
        return totalIntensity;
    }
    public void sortPoints(boolean sortByIntensity) {
        points = getSortedPoints(sortByIntensity);
    }

    /**
     * Return a sorted list (incremental by m2z or intensity)
     * for the Points in this PointList
     * @param sortByIntensity - indicate how the list should be sorted,
     *                          true for sort by intensity,
     *                          false for sort by M2z
     * @note user must not modify the List returned. 
     * 
     */
    public synchronized ArrayList<Point> getSortedPoints(boolean sortByIntensity) {
        ArrayList sortedPoints =  new ArrayList(DEFAULTNUMPEAKS);
        if (sortByIntensity) {
            if (pointsSortedByIntensity != null) {
                return pointsSortedByIntensity;
            } else {
                pointsSortedByIntensity = sortedPoints;
            }
        } else {
            if (pointsSortedByXValue != null) {
                return pointsSortedByXValue;
            } else {
                pointsSortedByXValue = sortedPoints;
            }
        }
        for (Point p : points) {
            sortedPoints.add(p);
        }

        Collections.sort(sortedPoints, new PointComparator(sortByIntensity));
        return sortedPoints; 
    }

    public ListIterator<Point> getPoints() {
        return points.listIterator();
    }

    public String getListType()
    {
        return listType;
    }

    public void setListType(String listType)
    {
        this.listType = listType;
    }

    public double getMaxXValue() {
        List<Point> sortedPoints = getSortedPoints(SORTBYXVALUE);
        return sortedPoints.get(numPoints()-1).getXValue();
    } 
    public double getMinXValue() {
        List<Point> sortedPoints = getSortedPoints(SORTBYXVALUE);
        return sortedPoints.get(0).getXValue();
    } 
    public double getMaxIntensity() {
        List<Point> sortedPoints = getSortedPoints(SORTBYINTENSITY);
        return sortedPoints.get(points.size()-1).getIntensity();
    }
    public double getMinIntensity() {
        List<Point> sortedPoints = getSortedPoints(SORTBYINTENSITY);
        return sortedPoints.get(0).getIntensity();
    }
}
