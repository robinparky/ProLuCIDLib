
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
public class ResidueParams {

    private String modificationMass;
    private String potentialModificationMass;
    private String potentialModificationMotif;
    
    public ResidueParams() {

    }

    public String getModificationMass() {
        return modificationMass;
    }

    public void setModificationMass(String modificationMass) {
        this.modificationMass = modificationMass;
    }

    public String getPotentialModificationMass() {
        return potentialModificationMass;
    }

    public void setPotentialModificationMass(String potentialModificationMass) {
        this.potentialModificationMass = potentialModificationMass;
    }

    public String getPotentialModificationMotif() {
        return potentialModificationMotif;
    }

    public void setPotentialModificationMotif(String potentialModificationMotif) {
        this.potentialModificationMotif = potentialModificationMotif;
    }
}
