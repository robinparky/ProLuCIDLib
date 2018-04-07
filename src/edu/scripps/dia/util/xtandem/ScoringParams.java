
/*
* Copyright (c) 2008 Integrated Proteomics Applications.  All rights reserved.  
*/

package edu.scripps.dia.util.xtandem;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email robinparky@yahoo.com
 * Created on Sep 1, 2009 
 * $Revision: 1.1.1.1 $
 * $Date: 2010/12/19 23:44:06 $
 */
public class ScoringParams {

    private String minimumIonCount = "4";
    private String maximumMissedCleavageSites="1";
    private String xIons="no";
    private String yIons="yes";
    private String zIons="no";
    private String aIons="no";
    private String bIons="yes";
    private String cIons="no";
    private String cyclicPermutation="no";
    private String includeReverse="no";
    
    public ScoringParams() {

    }

    public String getMinimumIonCount() {
        return minimumIonCount;
    }

    public void setMinimumIonCount(String minimumIonCount) {
        this.minimumIonCount = minimumIonCount;
    }

    public String getMaximumMissedCleavageSites() {
        return maximumMissedCleavageSites;
    }

    public void setMaximumMissedCleavageSites(String maximumMissedCleavageSites) {
        this.maximumMissedCleavageSites = maximumMissedCleavageSites;
    }

    public String getXIons() {
        return xIons;
    }

    public void setXIons(String xIons) {
        this.xIons = xIons;
    }

    public String getYIons() {
        return yIons;
    }

    public void setYIons(String yIons) {
        this.yIons = yIons;
    }

    public String getZIons() {
        return zIons;
    }

    public void setZIons(String zIons) {
        this.zIons = zIons;
    }

    public String getAIons() {
        return aIons;
    }

    public void setAIons(String aIons) {
        this.aIons = aIons;
    }

    public String getBIons() {
        return bIons;
    }

    public void setBIons(String bIons) {
        this.bIons = bIons;
    }

    public String getCIons() {
        return cIons;
    }

    public void setCIons(String cIons) {
        this.cIons = cIons;
    }

    public String getCyclicPermutation() {
        return cyclicPermutation;
    }

    public void setCyclicPermutation(String cyclicPermutation) {
        this.cyclicPermutation = cyclicPermutation;
    }

    public String getIncludeReverse() {
        return includeReverse;
    }

    public void setIncludeReverse(String includeReverse) {
        this.includeReverse = includeReverse;
    }
}
