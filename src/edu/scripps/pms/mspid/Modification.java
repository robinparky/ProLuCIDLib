/**
 * @file Modification.java
 * This is the source file for edu.scripps.pms.mspid.Modification
 * @author Tao Xu
 * @date $Date: 2007/03/09 20:54:16 $
 */

package edu.scripps.pms.mspid;


public class Modification {
    private char residue;
    private double massShift;
    private char symbol;
    private int numMods = 1;
    int hashcode;
//    private ArrayList<Double> differentialModifications = new ArrayList<Double>();

    public Modification(char residue, double massShift) {
        this.residue = residue;
        this.massShift = massShift;
        hashcode = new Double(residue + massShift).hashCode();
    }
    public Modification(char residue, double massShift, char symbol) {
        this.residue = residue;
        this.massShift = massShift;
        this.symbol = symbol;
        hashcode = new Double(residue + massShift).hashCode();
    }
    public double getMassShift() {
        return massShift;
    }
    public char getResidue() {
        return residue;
    }
    public void setSymbol(char s) {
        symbol = s;
    }
    public char getSymbol() {
        return symbol;
    }

    public boolean equals(Modification o) {
        return residue == o.residue && massShift == o.massShift;
    }
    public int hashCode() {
        return hashcode; 
    }
}



