
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

import java.io.File;

public class XTandemParams {

    
    private Integer proteinDbId;    
    private String databaseName = "";    
    private String databasePath = "";    
    
    private SpectrumParams spectrumParams;
    private OutputParams outputParams;
    private ProteinParams proteinParams;
    private RefineParams refineParams;
    private ResidueParams residueParams;
    private ScoringParams scoringParams;
    
    private String diffmods;
    
    
    private String addCTerminus="0";
    private String addNTerminus="0";
    private String addG="0";
    private String addA="0";
    private String addS="0";
    private String addP="0";
    private String addV="0";
    private String addT="0";
    private String addC="57.02146";
    private String addL="0";
    private String addI="0";
    private String addX="0";
    private String addN="0";
    private String addO="0";
    private String addB="0";
    private String addD="0";
    private String addQ="0";
    private String addK="0";
    private String addZ="0";
    private String addE="0";
    private String addM="0";
    private String addH="0";
    private String addF="0";
    private String addR="0";
    private String addY="0";
    private String addW="0";
    
    public XTandemParams() {

    }
    
    public String getEnzymeInfo() {
	StringBuffer sb = new StringBuffer();
	sb.append("[SEQUEST_ENZYME_INFO]").append("\n");
	sb.append("0.  No_Enzyme              0      -           -").append("\n");
	sb.append("1.  Trypsin                1      KR          P").append("\n");
	sb.append("2.  Chymotrypsin           1      FWY         P").append("\n");
	sb.append("3.  Clostripain            1      R           -").append("\n");
	sb.append("4.  Cyanogen_Bromide       1      M           -").append("\n");
	sb.append("5.  IodosoBenzoate         1      W           -").append("\n");
	sb.append("6.  Proline_Endopept       1      P           -").append("\n");
	sb.append("7.  Staph_Protease         1      E           -").append("\n");
	sb.append("8.  Trypsin_K              1      K           P").append("\n");
	sb.append("9.  Trypsin_R              1      R           P").append("\n");
	sb.append("10. AspN                   0      D           -").append("\n");
	sb.append("11. Cymotryp/Modified      1      FWYL        P").append("\n");
	sb.append("12. Elastase               1      ALIV        P").append("\n");
	sb.append("13. Elastase/Tryp/Chymo    1      ALIVKRWFY   P").append("\n");
	sb.append("").append("\n");

	return sb.toString();
    }
        
    public void readParams(String file) throws Exception {

	File f = new File(file);
	if(!f.exists())
	    return;
        
        /*

	Hashtable<String, String> ht = new Hashtable<String, String>();

	BufferedReader br = new BufferedReader(new FileReader(file));
	String eachLine = br.readLine();

	if(eachLine.contains("heavy_search_mass_shift")) {
	    String[] tmpArr = eachLine.split("=");
	    this.heavySearchMassShift = tmpArr[1].trim();
	    this.heavySearchMassShiftArr = heavySearchMassShift.split(" ");
	    this.labelSearchMass = true;
	}

	while ( (eachLine = br.readLine()) != null)// some H lines are followed by space, not tab. So, don't use H\t. :(
	{
	    eachLine = eachLine.trim();
	    if(eachLine.startsWith("#") || eachLine.startsWith("[") || "".equals(eachLine))
		continue;

	    String[] arr = eachLine.split("=");

	    if(arr.length<2)
		continue;

	    String[] tmpArr = arr[1].split(";");
	    ht.put(arr[0].trim(), tmpArr[0].trim());
	}

	if(null != ht.get("ion_cutoff_percentage"))
	    this.ionCutoffPercentage=ht.get("ion_cutoff_percentage").toString();
	if(null != ht.get("enzyme_number"))
	    this.enzymeNumber=ht.get("enzyme_number").toString();
	if(null != ht.get("add_C_Cysteine"))
	    this.addC=ht.get("add_C_Cysteine").toString();
	if(null != ht.get("num_description_lines"))
	    this.numDescriptionLines=ht.get("num_description_lines").toString();
	    this.massTypeParent=tmpValue;
	    if("1".equals(massTypeParent))
		this.resolution = "high";
	    else
		this.resolution = "low";
	}
	if(null != ht.get("protein_mass_filter"))
	    this.proteinMassFilter=ht.get("protein_mass_filter").toString();
	if(null != ht.get("diff_search_options"))
	    this.diffSearchOptions=ht.get("diff_search_options").toString();
	*/
    
    }

    public Integer getProteinDbId() {
        return proteinDbId;
    }

    public void setProteinDbId(Integer proteinDbId) {
        this.proteinDbId = proteinDbId;
    }

    public String getDatabaseName() {
        return databaseName;
    }

    public void setDatabaseName(String databaseName) {
        this.databaseName = databaseName;
    }

    public String getDatabasePath() {
        return databasePath;
    }

    public void setDatabasePath(String databasePath) {
        this.databasePath = databasePath;
    }

    public SpectrumParams getSpectrumParams() {
        return spectrumParams;
    }

    public void setSpectrumParams(SpectrumParams spectrumParams) {
        this.spectrumParams = spectrumParams;
    }

    public OutputParams getOutputParams() {
        return outputParams;
    }

    public void setOutputParams(OutputParams outputParams) {
        this.outputParams = outputParams;
    }

    public ProteinParams getProteinParams() {
        return proteinParams;
    }

    public void setProteinParams(ProteinParams proteinParams) {
        this.proteinParams = proteinParams;
    }

    public RefineParams getRefineParams() {
        return refineParams;
    }

    public void setRefineParams(RefineParams refineParams) {
        this.refineParams = refineParams;
    }

    public ResidueParams getResidueParams() {
        return residueParams;
    }

    public void setResidueParams(ResidueParams residueParams) {
        this.residueParams = residueParams;
    }

    public ScoringParams getScoringParams() {
        return scoringParams;
    }

    public void setScoringParams(ScoringParams scoringParams) {
        this.scoringParams = scoringParams;
    }

    public String getDiffmods() {
        return diffmods;
    }

    public void setDiffmods(String diffmods) {
        this.diffmods = diffmods;
    }

    public String getAddCTerminus() {
        return addCTerminus;
    }

    public void setAddCTerminus(String addCTerminus) {
        this.addCTerminus = addCTerminus;
    }

    public String getAddNTerminus() {
        return addNTerminus;
    }

    public void setAddNTerminus(String addNTerminus) {
        this.addNTerminus = addNTerminus;
    }

    public String getAddG() {
        return addG;
    }

    public void setAddG(String addG) {
        this.addG = addG;
    }

    public String getAddA() {
        return addA;
    }

    public void setAddA(String addA) {
        this.addA = addA;
    }

    public String getAddS() {
        return addS;
    }

    public void setAddS(String addS) {
        this.addS = addS;
    }

    public String getAddP() {
        return addP;
    }

    public void setAddP(String addP) {
        this.addP = addP;
    }

    public String getAddV() {
        return addV;
    }

    public void setAddV(String addV) {
        this.addV = addV;
    }

    public String getAddT() {
        return addT;
    }

    public void setAddT(String addT) {
        this.addT = addT;
    }

    public String getAddC() {
        return addC;
    }

    public void setAddC(String addC) {
        this.addC = addC;
    }

    public String getAddL() {
        return addL;
    }

    public void setAddL(String addL) {
        this.addL = addL;
    }

    public String getAddI() {
        return addI;
    }

    public void setAddI(String addI) {
        this.addI = addI;
    }

    public String getAddX() {
        return addX;
    }

    public void setAddX(String addX) {
        this.addX = addX;
    }

    public String getAddN() {
        return addN;
    }

    public void setAddN(String addN) {
        this.addN = addN;
    }

    public String getAddO() {
        return addO;
    }

    public void setAddO(String addO) {
        this.addO = addO;
    }

    public String getAddB() {
        return addB;
    }

    public void setAddB(String addB) {
        this.addB = addB;
    }

    public String getAddD() {
        return addD;
    }

    public void setAddD(String addD) {
        this.addD = addD;
    }

    public String getAddQ() {
        return addQ;
    }

    public void setAddQ(String addQ) {
        this.addQ = addQ;
    }

    public String getAddK() {
        return addK;
    }

    public void setAddK(String addK) {
        this.addK = addK;
    }

    public String getAddZ() {
        return addZ;
    }

    public void setAddZ(String addZ) {
        this.addZ = addZ;
    }

    public String getAddE() {
        return addE;
    }

    public void setAddE(String addE) {
        this.addE = addE;
    }

    public String getAddM() {
        return addM;
    }

    public void setAddM(String addM) {
        this.addM = addM;
    }

    public String getAddH() {
        return addH;
    }

    public void setAddH(String addH) {
        this.addH = addH;
    }

    public String getAddF() {
        return addF;
    }

    public void setAddF(String addF) {
        this.addF = addF;
    }

    public String getAddR() {
        return addR;
    }

    public void setAddR(String addR) {
        this.addR = addR;
    }

    public String getAddY() {
        return addY;
    }

    public void setAddY(String addY) {
        this.addY = addY;
    }

    public String getAddW() {
        return addW;
    }

    public void setAddW(String addW) {
        this.addW = addW;
    }

}
