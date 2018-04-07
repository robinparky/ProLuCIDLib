
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
public class RefineParams {

    private String modificationMass;
    private String sequencePath;
    private String ticPercent = "20";
    private String spectrumSynthesis = "yes";
    private String maximumValidExpectationValue = "0.1";
    private String potentialNterminusModifications = "+42.010565@[";
    private String potentialCterminusModifications;
    private String unanticipatedCleavage = "yes";
    private String potentialModificationMass;
    private String pointMutations = "no";
    private String usePotentialModificationsForFullRefinement = "no";
    //private String pointMutations ;
    private String potentialModificationMotif;

    public RefineParams() {

    }

    public String getModificationMass() {
        return modificationMass;
    }

    public void setModificationMass(String modificationMass) {
        this.modificationMass = modificationMass;
    }

    public String getSequencePath() {
        return sequencePath;
    }

    public void setSequencePath(String sequencePath) {
        this.sequencePath = sequencePath;
    }

    public String getTicPercent() {
        return ticPercent;
    }

    public void setTicPercent(String ticPercent) {
        this.ticPercent = ticPercent;
    }

    public String getSpectrumSynthesis() {
        return spectrumSynthesis;
    }

    public void setSpectrumSynthesis(String spectrumSynthesis) {
        this.spectrumSynthesis = spectrumSynthesis;
    }

    public String getMaximumValidExpectationValue() {
        return maximumValidExpectationValue;
    }

    public void setMaximumValidExpectationValue(String maximumValidExpectationValue) {
        this.maximumValidExpectationValue = maximumValidExpectationValue;
    }

    public String getPotentialNterminusModifications() {
        return potentialNterminusModifications;
    }

    public void setPotentialNterminusModifications(String potentialNterminusModifications) {
        this.potentialNterminusModifications = potentialNterminusModifications;
    }

    public String getPotentialCterminusModifications() {
        return potentialCterminusModifications;
    }

    public void setPotentialCterminusModifications(String potentialCterminusModifications) {
        this.potentialCterminusModifications = potentialCterminusModifications;
    }

    public String getUnanticipatedCleavage() {
        return unanticipatedCleavage;
    }

    public void setUnanticipatedCleavage(String unanticipatedCleavage) {
        this.unanticipatedCleavage = unanticipatedCleavage;
    }

    public String getPotentialModificationMass() {
        return potentialModificationMass;
    }

    public void setPotentialModificationMass(String potentialModificationMass) {
        this.potentialModificationMass = potentialModificationMass;
    }

    public String getPointMutations() {
        return pointMutations;
    }

    public void setPointMutations(String pointMutations) {
        this.pointMutations = pointMutations;
    }

    public String getUsePotentialModificationsForFullRefinement() {
        return usePotentialModificationsForFullRefinement;
    }

    public void setUsePotentialModificationsForFullRefinement(String usePotentialModificationsForFullRefinement) {
        this.usePotentialModificationsForFullRefinement = usePotentialModificationsForFullRefinement;
    }

    public String getPotentialModificationMotif() {
        return potentialModificationMotif;
    }

    public void setPotentialModificationMotif(String potentialModificationMotif) {
        this.potentialModificationMotif = potentialModificationMotif;
    }
}
