/**
 * @file DiffMod.java
 * This is the source file for edu.scripps.pms.mspid.Modification
 * @author Tao Xu
 * @date $Date: 2009/07/05 05:34:34 $
 */

package edu.scripps.pms.mspid;


import java.text.DecimalFormat;


public class DiffMod {
    private double massShift;
    private char symbol;
    //private HashSet<Modification> mods = new HashSet<Modification>();
    int hashcode;
    private String info;
    private String modifiableResidues = null;
    private boolean [] modifiables = new boolean[256];
    private double dbinwidthMassShift;
    private double neutralLoss;
    private double dbinwidthNeutralLoss;
    private int intMassShift;    
    public static DecimalFormat threeDigits = new DecimalFormat("0.000");
//    private ArrayList<Double> differentialModifications = new ArrayList<Double>();
    public String getModifiableResidues() {
        if(modifiableResidues == null) {
            StringBuffer sb = new StringBuffer(20);
            for(int i = 0; i < MassSpecConstants.NUMCHARS; i++) {
                if(modifiables[i]) {
                    sb.append((char)i);
                }
            }
            modifiableResidues = sb.toString();
        }

        return modifiableResidues;
        
    }

    public String toString() {
        StringBuffer sb = new StringBuffer(100);
        //sb.append(info);
        //sb.append("\t");
        for(int i = 0; i < MassSpecConstants.NUMCHARS; i++) {
            if(modifiables[i]) {
                sb.append((char)i);
            }
        }
        sb.append(symbol);
        sb.append("=");
        if(massShift > 0) {
            sb.append("+");
        } 
        sb.append(massShift + " ");
        return sb.toString();
    }
    public DiffMod(double massShift, char symbol) {
        this.symbol = symbol;
        this.massShift = massShift;
        intMassShift = (int)(massShift* MassCalculator.MASSACCURACYFACTOR+0.5);
        dbinwidthMassShift = massShift* MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        info = "(" + threeDigits.format(massShift) + ")";
    }
    public DiffMod(double massShift, double neutralLoss, char symbol) {
        this.symbol = symbol;
        this.massShift = massShift;
    
        dbinwidthMassShift = massShift* MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        this.neutralLoss = neutralLoss;  
        dbinwidthNeutralLoss = neutralLoss* MassSpecConstants.DBINWIDTH;
        info = "(" + threeDigits.format(massShift) + ")";
    }
    public void setModifiable(byte residue, boolean modifiable) {
        modifiables[residue] = modifiable;
    }
    public void setModifiable(boolean [] residues) {
        modifiables = residues;
    }
    public boolean isModifiable(byte residue) {
        return modifiables[residue];
    }
    
    public double getDbinwidthNeutralLoss() {
        return dbinwidthNeutralLoss;
    }
    public double getDbinwidthMassShift() {
        return dbinwidthMassShift;
    }
    public void setMassShift(double shift) {
        massShift = shift;
        intMassShift = (int)(massShift* MassCalculator.MASSACCURACYFACTOR+0.5);
        dbinwidthMassShift = massShift* MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        info = "(" + threeDigits.format(massShift) + ")";
        
    }
    public double getMassShift() {
        return massShift;
    }
    public double getNeutralLoss() {
        return neutralLoss;
    }
    public int getIntMassShift() {
        return intMassShift;
    }
    public char getSymbol() {
        return symbol;
    }
    public String getModInfo() {
        return info;
        //return ""+symbol;
    }
    public void setModInfo(String s) {
        info = s;
    }
    public boolean equals(DiffMod o) {
        return massShift == o.massShift;
    }
    public int hashCode() {
        return hashcode; 
    }
}



