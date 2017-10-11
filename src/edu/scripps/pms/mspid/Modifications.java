/**
 * @file Modification.java
 * This is the source file for edu.scripps.pms.mspid.Modification
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:41:49 $
 */

package edu.scripps.pms.mspid;


import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;


public class Modifications {
    private double staticCTermMod;
    private double staticNTermMod;
    // diff mods
    private TerminalModification nTerm;  // nterm diff mod
    private TerminalModification cTerm;
    private HashSet<Modification> staticMods = new HashSet<Modification>(20);
    private ArrayList<DiffMod> diffMods = new ArrayList<DiffMod>(20);
    private HashSet<DiffMod>[] residue2DiffMods = (HashSet<DiffMod>[]) new HashSet[256];
    private boolean modifiables [] = new boolean[256];
    private String info = null;
    
    public Modifications copy() {
        Modifications m = new Modifications();
        m.staticCTermMod = staticCTermMod;
        m.staticNTermMod = staticNTermMod;
        m.nTerm = nTerm;
        m.cTerm = cTerm;
        // static mods are ignored
        for(Iterator<DiffMod> it = getDiffMods(); it.hasNext();) {
            m.addDiffMod(it.next());
        }
        return m;
    }

    // num of terminal diff mods does not count
    public int getNumDiffMods() {
        /*
        int numDiffMods = 0;
        if(cTerm != null) {
            numDiffMods++;
        }
        if(nTerm != null) {
            numDiffMods++;
        }
        return numDiffMods + diffMods.size();
        */
        return diffMods.size();
    }
    private void init() {

        for(int i = 0; i < residue2DiffMods.length; i++) {
            residue2DiffMods[i] = new HashSet<DiffMod>();
        }
    }
    public Modifications() {
        init();
    }
    public Modifications(TerminalModification t, boolean isNTerm) {
        init();
        if(isNTerm) {
            setNTerm(t);
        } else {
            setCTerm(t);
        }
    }
    
    public Modifications(TerminalModification n, TerminalModification c) {
        init();
        setNTerm(n);
        setCTerm(c);
    }
    public void setNTerm(TerminalModification tm) {
        nTerm = tm;
    }
    public void setCTerm(TerminalModification tm) {
        cTerm = tm;
    }
    public void setStaticNTermMod(double m) {
        staticNTermMod = m;
    }
    public void setStaticCTermMod(double m) {
        staticCTermMod = m;
    }
    public double getStaticTerminalMods() {
        return staticNTermMod + staticCTermMod;
    }    
    public double getMassShift() {
        double massShift = 0;
        for(Iterator<DiffMod> it = diffMods.iterator(); it.hasNext();) {
            DiffMod d = it.next();
            massShift += d.getMassShift();
        }
        return /*getNTermMassShift() + getCTermMassShift() +*/ massShift;
    }
    public double getInternalDiffModsMass() {
        double result = 0; 
        for(Iterator<DiffMod> it = diffMods.iterator(); it.hasNext();) {
            result += it.next().getMassShift();
        }
        return result;
    }
    public double getDiffModsShift() {
        //double result = getStaticTerminalMods();
        double result = 0;
        result += getDiffNTermMod();
        result += getDiffCTermMod();
        for(Iterator<DiffMod> it = diffMods.iterator(); it.hasNext();) {
            result += it.next().getMassShift();
        }
        return result;
    }
    public void addDiffMod(DiffMod m) {
//System.out.println("in addDiffMod massShift: " + m.getMassShift());
       if(m != null && m.getMassShift() != 0) {
           diffMods.add(m);
           for(byte c = 0; c < 127; c++) {
               if(m.isModifiable(c)) {
                   residue2DiffMods[c].add(m);
                   modifiables[c] = true;
               }
           } 
       }
       
    }
    public HashSet<DiffMod> residue2DiffMods(byte residue) {
        return residue2DiffMods[residue];
    }
    public boolean isModifiable(byte residue) {
//System.out.println("in ismodifiable residue is: " + residue);

        return modifiables[residue];
    }
    public void addStaticMod(Modification m) {
        staticMods.add(m);
    }
    public double getDiffNTermMod() {
        return nTerm == null? 0 : nTerm.getMassShift();
    }
    public double getDiffCTermMod() {
        return cTerm == null? 0 : cTerm.getMassShift();
    }
    public double getNTermMassShift() {
        return nTerm == null? staticNTermMod : nTerm.getMassShift() + staticNTermMod;
    }
    public double getCTermMassShift() {
        return cTerm == null? staticCTermMod : cTerm.getMassShift() + staticCTermMod;
    }
    public double getStaticNTermMod() {
        return staticNTermMod;
    }
   
    public double getStaticCTermMod() {
        return staticCTermMod;
    }
    public TerminalModification getCTerm() {
        return cTerm;
    }
    public TerminalModification getNTerm() {
        return nTerm;
    }
    public Iterator<DiffMod> getDiffMods() {
        return diffMods.iterator();
    }
    public String getInfo() {
        if(info == null) {        
            StringBuffer sb = new StringBuffer(200);
            sb.append("cTerm: " + getCTermMassShift() + "\tnTerm: " + getNTermMassShift() + "\tnumDiffMods: " + diffMods.size());
            for(Iterator<DiffMod> it = diffMods.iterator(); it.hasNext();) {
                sb.append("\t" + it.next().getModInfo());
            }
            info = sb.toString();
        }
        return info;
    }
    public Iterator<Modification> getStaticMods() {
        return staticMods.iterator();
    }

    public boolean equals(Object o) {
        if(o == null) {
            return false;
        }
        Modifications m = (Modifications) o;
//System.out.println(getInfo() + "\tm info: " + m.getInfo());
        boolean allTheSame = true;
        allTheSame = (getNTermMassShift() == m.getNTermMassShift() &&
               getCTermMassShift() == m.getCTermMassShift() &&
               getMassShift() == m.getMassShift() &&
               getNumDiffMods() == m.getNumDiffMods());                
        if(!allTheSame) {
            return false;
        }  
        HashSet<DiffMod> mySet = new HashSet<DiffMod>(diffMods);
        HashSet<DiffMod> oSet = new HashSet<DiffMod>(m.diffMods);
        return mySet.equals(oSet);
        
    }
    public int hashCode() {
        return getInfo().hashCode();
    }
}



