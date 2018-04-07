
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
public class OutputParams {

    private String logPath;
    private String message;
    private String oneSequenceCopy="no";
    private String sequencePath;
    private String path;
    private String sortResultsBy;
    private String pathHashing;
    private String xslPath;
    private String parameters;
    private String performance;
    private String spectra;
    private String histograms;
    private String proteins;
    private String sequences;
    //private String oneSequenceCopy;
    private String results;
    private String maximumValidExpectationValue;
    private String histogramColumnWidth;
    
    public OutputParams() {

    }

    public String getLogPath() {
        return logPath;
    }

    public void setLogPath(String logPath) {
        this.logPath = logPath;
    }

    public String getMessage() {
        return message;
    }

    public void setMessage(String message) {
        this.message = message;
    }

    public String getOneSequenceCopy() {
        return oneSequenceCopy;
    }

    public void setOneSequenceCopy(String oneSequenceCopy) {
        this.oneSequenceCopy = oneSequenceCopy;
    }

    public String getSequencePath() {
        return sequencePath;
    }

    public void setSequencePath(String sequencePath) {
        this.sequencePath = sequencePath;
    }

    public String getPath() {
        return path;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getSortResultsBy() {
        return sortResultsBy;
    }

    public void setSortResultsBy(String sortResultsBy) {
        this.sortResultsBy = sortResultsBy;
    }

    public String getPathHashing() {
        return pathHashing;
    }

    public void setPathHashing(String pathHashing) {
        this.pathHashing = pathHashing;
    }

    public String getXslPath() {
        return xslPath;
    }

    public void setXslPath(String xslPath) {
        this.xslPath = xslPath;
    }

    public String getParameters() {
        return parameters;
    }

    public void setParameters(String parameters) {
        this.parameters = parameters;
    }

    public String getPerformance() {
        return performance;
    }

    public void setPerformance(String performance) {
        this.performance = performance;
    }

    public String getSpectra() {
        return spectra;
    }

    public void setSpectra(String spectra) {
        this.spectra = spectra;
    }

    public String getHistograms() {
        return histograms;
    }

    public void setHistograms(String histograms) {
        this.histograms = histograms;
    }

    public String getProteins() {
        return proteins;
    }

    public void setProteins(String proteins) {
        this.proteins = proteins;
    }

    public String getSequences() {
        return sequences;
    }

    public void setSequences(String sequences) {
        this.sequences = sequences;
    }

    public String getResults() {
        return results;
    }

    public void setResults(String results) {
        this.results = results;
    }

    public String getMaximumValidExpectationValue() {
        return maximumValidExpectationValue;
    }

    public void setMaximumValidExpectationValue(String maximumValidExpectationValue) {
        this.maximumValidExpectationValue = maximumValidExpectationValue;
    }

    public String getHistogramColumnWidth() {
        return histogramColumnWidth;
    }

    public void setHistogramColumnWidth(String histogramColumnWidth) {
        this.histogramColumnWidth = histogramColumnWidth;
    }
}
