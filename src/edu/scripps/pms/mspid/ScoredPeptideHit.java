/**
 * @file ScoredPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.ScoredPeptideHit
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:41:49 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.dia.LibrarySpectra;

import java.util.LinkedList;
import java.util.List;

public class ScoredPeptideHit implements Comparable<ScoredPeptideHit> {
    // Each ScoredPeptideHit may be associated with multiple PeptideHit
    private double xcorr;
    private double zscore;
    private double pScore; // probability score
    private double pValue; // p-value based on xcorr
    private double eValue; // e-value based on xcorr
    private int xcorrRank;
    private int pScoreRank;
    private LinkedList<PeptideHit> hits = new LinkedList<PeptideHit>();
    private String sequence;
    private double prcMass;
    private double rSqaured;
    private double retTime;
    private int numPeaks=-1;
    private int numMatchedPeaks = -1;
    private boolean isDecoy = false;
    private int scanHi;
    private int scanLow;
    private String ms2Filename;
    
    // 0 for sort by binomial probability,
    // 1 for sort by xcorr 
    private int primaryScoreType;
    private int secondaryScoreType;
    private boolean hasPeptide = true;
    public void setZscore(double z) {
        zscore = z;
    }
    public double getZscore() {
        return zscore;
    }
    public int [] getTheorMasses() {
        return hits.get(0).getTheorMasses();
    }
    public double getTheorMass() {
        return hits.get(0).getTheorMass();
    }
    private LibrarySpectra librarySpectra;

    public ScoredPeptideHit(double probability) {
        hasPeptide = false;
        pScore = probability;
        primaryScoreType = 0;
    }

    public ScoredPeptideHit(String seq, int primaryscoretype, int secondaryscoretype) {
        sequence = seq;
        this.primaryScoreType = primaryscoretype;
        this.secondaryScoreType = secondaryscoretype;
    }
    public void setPValue(double p) {
        pValue = p;
    }
    public double setPValue() {
        return pValue;
    }
    public void setEValue(double e) {
        eValue = e;
    }
    public double eetEValue() {
        return eValue;
    }
    public boolean isModified() {
        return hits.get(0).isModified();
    }
    public String getExactSequence() {
        return hits.get(0).getExactSequence();
    }
    // get AA sequence without modification symbols
    public String getOriginalSequence() {
        
        return hits.get(0).getSequence();
    }
    public String getSequence() {
        return librarySpectra.sequence;
        //return hits.get(0).getExtendedSequence();
    }
    public String getExtendedSequence() {
        return librarySpectra.sequence;
    }
    public int getXCorrRank() {
        return xcorrRank;
    }
    public int getPScoreRank() {
        return pScoreRank;
    }
    public int getSecondaryRank() {
        switch(primaryScoreType) {
            case 0: return xcorrRank;
            case 1: return pScoreRank; 
            case 2: return pScoreRank; 
            default: return xcorrRank;
        }    
    }
    public double getSecondaryScore() {
        switch(secondaryScoreType) {
            case 0: return pScore;
            case 1: return xcorr; 
            case 2: return zscore; 
            case 3: return pValue; 
            case 4: return eValue; 
            default: return xcorr;
        }    
    }
    public int getPrimaryRank() {
        switch(primaryScoreType) {
            case 0: return pScoreRank;
            case 1: return xcorrRank; 
            case 2: return xcorrRank; 
            default: return pScoreRank;
        }    
    }
    public double getPrimaryScore() {
        return pScore;

     /*   switch(primaryScoreType) {
            case 0: return pScore;
            case 1: return xcorr; 
            case 2: return zscore; 
            default: return pScore;
        }    */
    }
    public void setPrimaryRank(int r) {
        
        switch(primaryScoreType) {
            case 0: pScoreRank = r;
            case 1: xcorrRank = r; 
            default: pScoreRank = r;
        }    
    }
    public void setXCorr(double x) {
        xcorr = x;
    }
    public void setPScore(double x) {
        pScore = x;
    }
    public void setPScoreRank(int r) {
        pScoreRank = r;
    }
    public void setXCorrRank(int r) {
        xcorrRank = r;
    }
    public double getPScore() {
        return pScore;
    }
    public double getXCorr() {
        return xcorr;
    }
    public int getNumPeaks() {
        return numPeaks;
        //return hits.get(0).getNumPeaks();
    }
    public int getNumPeaksMatched() {
       return numMatchedPeaks;
       // return hits.get(0).getNumPeaksMatched();
    } 
    public void addPeptideHit(PeptideHit p) {
        hits.add(p);
    } 
    public List<PeptideHit> getPeptideHits() {
        return hits;
    }

    // sort by xcorr 
    public int compareTo(ScoredPeptideHit s) {
        double f = s.xcorr- this.xcorr;  // difference

        if (f > 0) {
            return 1;
        } else if (f < 0) {
            return -1;
        } else {
            if(s.isModified() == this.isModified()) {
                return 0;
            } else if(s.isModified()) {
                return -1;
            } else {
                return 1;
            }
        }

    }
    public DiffMod[] getDiffMods()
    {
        return hits.get(0).getDiffMods();
    }
    public Modifications getModifications()
    {
        return hits.get(0).getModifications();
    }

    public boolean hasPeptide() {
        return hasPeptide;
    }

    public void setMs2CompareValues(int hiScan, int lowScan, String filename) {
        ms2Filename = filename;
        scanHi = hiScan;
        scanLow = lowScan;
    }

    public int getScanHi() {
        return scanHi;
    }

    public int getScanLow() {
        return scanLow;
    }

    public String getMs2Filename() {
        return ms2Filename;
    }

    public void setIsDecoy(boolean isDecoy) {
        this.isDecoy = isDecoy;
    }
    public boolean isDecoy()
    {
        return  isDecoy;
    }

    public double getPrcMass()
    {
        return prcMass;
    }

    public void setPrcMass(double prcMass)
    {
        this.prcMass = prcMass;
    }

    public double getrSqaured() {
        return rSqaured;
    }

    public void setrSqaured(double rSqaured) {
        this.rSqaured = rSqaured;
    }

    public double getRetTime() {
        return retTime;
    }

    public void setRetTime(double retTime) {
        this.retTime = retTime;
    }

    public void setNumPeaks(int numPeaks) {
        this.numPeaks = numPeaks;
    }

    public void setNumMatchedPeaks(int numMatchedPeaks) {
        this.numMatchedPeaks = numMatchedPeaks;
    }

    public LibrarySpectra getLibrarySpectra() {
        return librarySpectra;
    }

    public void setLibrarySpectra(LibrarySpectra librarySpectra) {
        this.librarySpectra = librarySpectra;
    }
}
