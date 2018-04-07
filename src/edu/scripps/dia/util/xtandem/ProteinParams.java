
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
public class ProteinParams {

    private String taxon="default";
    private String cleavageSite;
    private String modifiedResidueMassFile;
    private String cleavageCterminalMassChange;
    private String cleavageNterminalMassChange;
    private String nterminalResidueModificationMass;
    private String cterminalResidueModificationMass;
    private String homologManagement = "no";
    private Integer enzymeNumber;
    
    
    public ProteinParams() {

    }

    public String getTaxon() {
        return taxon;
    }

    public void setTaxon(String taxon) {
        this.taxon = taxon;
    }

    public String getCleavageSite() {
        return cleavageSite;
    }
    
    public void setCleavageSite(String cleavageSite) {
        this.cleavageSite = cleavageSite;
    }

    public String getModifiedResidueMassFile() {
        return modifiedResidueMassFile;
    }

    public void setModifiedResidueMassFile(String modifiedResidueMassFile) {
        this.modifiedResidueMassFile = modifiedResidueMassFile;
    }

    public String getCleavageCterminalMassChange() {
        return cleavageCterminalMassChange;
    }

    public void setCleavageCterminalMassChange(String cleavageCterminalMassChange) {
        this.cleavageCterminalMassChange = cleavageCterminalMassChange;
    }

    public String getCleavageNterminalMassChange() {
        return cleavageNterminalMassChange;
    }

    public void setCleavageNterminalMassChange(String cleavageNterminalMassChange) {
        this.cleavageNterminalMassChange = cleavageNterminalMassChange;
    }

    public String getNterminalResidueModificationMass() {
        return nterminalResidueModificationMass;
    }

    public void setNterminalResidueModificationMass(String nterminalResidueModificationMass) {
        this.nterminalResidueModificationMass = nterminalResidueModificationMass;
    }

    public String getCterminalResidueModificationMass() {
        return cterminalResidueModificationMass;
    }

    public void setCterminalResidueModificationMass(String cterminalResidueModificationMass) {
        this.cterminalResidueModificationMass = cterminalResidueModificationMass;
    }

    public String getHomologManagement() {
        return homologManagement;
    }

    public void setHomologManagement(String homologManagement) {
        this.homologManagement = homologManagement;
    }

    public Integer getEnzymeNumber() {
        return enzymeNumber;
    }

    public void setEnzymeNumber(Integer enzymeNumber) {
        this.enzymeNumber = enzymeNumber;
    }
    

    public String getEnzyme() {
        
        String enzyme=null;
        
        switch(this.enzymeNumber) {
            case 1:                 
                enzyme = "[RK]|{P}";
                break;
                
            case 2:
                enzyme = "[RK]|{}";
                break;
            
            case 3:
                enzyme = "[FWY]|{}";
                break;
                
            case 4:
                enzyme = "[R]|{}";
                break;
                
            case 5:
                enzyme = "[M]|{}";
                break;

            case 6:
                enzyme = "[W]|{}";
                break;

            case 7:
                enzyme = "[P]|{}";
                break;
            
            case 8:
                enzyme = "[E]|{}";
                break;
        
            case 9:
                enzyme = "[K]|{}";
                break;

            case 10:
                enzyme = "[R]|{}";
                break;
            case 11:
                enzyme = "[D]|{}";
                break;
            
            case 12:
                enzyme = "[FWYL]|{}";
                break;
        
            case 13:
                enzyme = "[ALIV]|{}";
                break;

            case 14:
                enzyme = "[ALIVKRWFY]|{}";
                break;
                
            default:
                enzyme = "";
                break;
                
        }
        
        return enzyme;

    }    

}
