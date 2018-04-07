
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
public class SpectrumParams {

    private String fragmentMonoisotopicMassError="0.4";
    private String parentMonoisotopicMassErrorPlus="100";
    private String parentMonoisotopicMassErrorMinus="100";
    private String parentMonoisotopicMassIsotopeError="yes";
    private String fragmentMonoisotopicMassErrorUnits="Daltons";
    private String parentMonoisotopicMassErrorUnits="ppm";
    private String fragmentMassType;
    private String dynamicRange = "100.0";
    private String totalPeaks = "50";
    private String maximumParentCharge = "4";
    private String useNoiseSuppression = "yes";
    private String minimumParentMh = "500.0";
    private String minimumFragmentMz = "150.0";
    private String minimumPeaks = "15";
    private String threads = "1";
    private String sequenceBatchSize = "1000";
   
    public SpectrumParams() {

    }

    public String getFragmentMonoisotopicMassError() {
        return fragmentMonoisotopicMassError;
    }

    public void setFragmentMonoisotopicMassError(String fragmentMonoisotopicMassError) {
        this.fragmentMonoisotopicMassError = fragmentMonoisotopicMassError;
    }

    public String getParentMonoisotopicMassErrorPlus() {
        
        return parentMonoisotopicMassErrorPlus;
    }

    public void setParentMonoisotopicMassErrorPlus(String parentMonoisotopicMassErrorPlus) {
        this.parentMonoisotopicMassErrorPlus = parentMonoisotopicMassErrorPlus;
    }

    public String getParentMonoisotopicMassErrorMinus() {
        return parentMonoisotopicMassErrorMinus;
    }

    public void setParentMonoisotopicMassErrorMinus(String parentMonoisotopicMassErrorMinus) {
        this.parentMonoisotopicMassErrorMinus = parentMonoisotopicMassErrorMinus;
    }

    public String getParentMonoisotopicMassIsotopeError() {
        return parentMonoisotopicMassIsotopeError;
    }

    public void setParentMonoisotopicMassIsotopeError(String parentMonoisotopicMassIsotopeError) {
        this.parentMonoisotopicMassIsotopeError = parentMonoisotopicMassIsotopeError;
    }

    public String getFragmentMonoisotopicMassErrorUnits() {
        return fragmentMonoisotopicMassErrorUnits;
    }

    public void setFragmentMonoisotopicMassErrorUnits(String fragmentMonoisotopicMassErrorUnits) {
        this.fragmentMonoisotopicMassErrorUnits = fragmentMonoisotopicMassErrorUnits;
    }

    public String getParentMonoisotopicMassErrorUnits() {
        return parentMonoisotopicMassErrorUnits;
    }

    public void setParentMonoisotopicMassErrorUnits(String parentMonoisotopicMassErrorUnits) {
        this.parentMonoisotopicMassErrorUnits = parentMonoisotopicMassErrorUnits;
    }

    public String getFragmentMassType() {
        return fragmentMassType;
    }

    public void setFragmentMassType(String fragmentMassType) {
        this.fragmentMassType = fragmentMassType;
    }

    public String getDynamicRange() {
        return dynamicRange;
    }

    public void setDynamicRange(String dynamicRange) {
        this.dynamicRange = dynamicRange;
    }

    public String getTotalPeaks() {
        return totalPeaks;
    }

    public void setTotalPeaks(String totalPeaks) {
        this.totalPeaks = totalPeaks;
    }

    public String getMaximumParentCharge() {
        return maximumParentCharge;
    }

    public void setMaximumParentCharge(String maximumParentCharge) {
        this.maximumParentCharge = maximumParentCharge;
    }

    public String getUseNoiseSuppression() {
        return useNoiseSuppression;
    }

    public void setUseNoiseSuppression(String useNoiseSuppression) {
        this.useNoiseSuppression = useNoiseSuppression;
    }

    public String getMinimumParentMh() {
        return minimumParentMh;
    }

    public void setMinimumParentMh(String minimumParentMh) {
        this.minimumParentMh = minimumParentMh;
    }

    public String getMinimumFragmentMz() {
        return minimumFragmentMz;
    }

    public void setMinimumFragmentMz(String minimumFragmentMz) {
        this.minimumFragmentMz = minimumFragmentMz;
    }

    public String getMinimumPeaks() {
        return minimumPeaks;
    }

    public void setMinimumPeaks(String minimumPeaks) {
        this.minimumPeaks = minimumPeaks;
    }

    public String getThreads() {
        return threads;
    }

    public void setThreads(String threads) {
        this.threads = threads;
    }

    public String getSequenceBatchSize() {
        return sequenceBatchSize;
    }

    public void setSequenceBatchSize(String sequenceBatchSize) {
        this.sequenceBatchSize = sequenceBatchSize;
    }
}
