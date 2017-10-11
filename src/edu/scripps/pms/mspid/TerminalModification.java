/**
 * @file TerminalModification.java
 * This is the source file for edu.scripps.pms.mspid.Modification
 * @author Tao Xu
 * @date $Date: 2007/12/02 18:37:12 $
 */

package edu.scripps.pms.mspid;


public class TerminalModification {
    private double massShift;
    private char symbol;
    private String info;
 
    public TerminalModification(char symbol, double massShift) {
        this.symbol = symbol;
        this.massShift = massShift;
        info = "(" + massShift + ")";
    }
    public double getMassShift() {
        return massShift;
    }
    public char getSymbol() {
        return symbol;
    }
    public String getModInfo() {
        return info;
    }
    public void setModInfo(String s) {
        info = s;
    }
    public String getInfoAsSymbol() {
        return "" + symbol;
    }

    public String getInfoAsMassShift() {
        return "(" + massShift + ")";
    }
}



