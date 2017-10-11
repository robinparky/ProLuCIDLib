/**
 * @file DeChargedModifiedPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.DeChargedModifiedPeptideHit
 * @author Tao Xu
 * @date $Date: 2009/12/07 18:50:11 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;

import java.util.ArrayList;
import java.util.Iterator;

public class DeChargedModifiedPeptideHit extends ModifiedPeptideHit {
     


    public DeChargedModifiedPeptideHit(Fasta parent, int start, int end, Modifications m) {
        super(parent, start, end, m);
    }

    public DeChargedModifiedPeptideHit(Peptide p, Modifications m) {
        super(p, m);

    }
    
    public DeChargedModifiedPeptideHit(ModifiedPeptideHit mph) {
        super(mph);

    }
     

    public DeChargedModifiedPeptideHit(Peptide p) {
        super(p);

    }



    public Iterator<ModifiedPeptideHit> getAllModifiedPeptideHits() {
        ArrayList<ModifiedPeptideHit> candidates = generateCandidates();
        ArrayList<ModifiedPeptideHit> dechargedCandidates = new ArrayList<ModifiedPeptideHit>();
        for(Iterator<ModifiedPeptideHit> it = candidates.iterator(); it.hasNext();) {
            ModifiedPeptideHit mph = it.next();
            dechargedCandidates.add((ModifiedPeptideHit)new DeChargedModifiedPeptideHit(mph));
        } 
        return dechargedCandidates.iterator();
    }
    public double getTheorMass() {
       
        double m = ppl.getMassCalculator().getPrecursorMass(peptide.getSequence());
        //m += internalDiffModMassAdded == 0? mods.getMassShift() : internalDiffModMassAdded;
        m += mods.getMassShift();
          
        //return ppl.getMassCalculator().getPrecursorMass(peptide.getSequence()) +  mods.getMassShift();
        return m;
    }
    // override the function in super class
    public double getNTermStart() { // b == n term
        if(mods != null) {
            return ppl.getNTermStart() + mods.getDiffNTermMod();
        } else {
            return ppl.getNTermStart();

        }
    }    
    // override the function in super class
    public double getNTermDbinwidthStart() { // b == n term
        if(mods != null) {
            return (ppl.getNTermStart() + mods.getDiffNTermMod())* MassSpecConstants.DBINWIDTH;
        } else {
            return ppl.getNTermStart()* MassSpecConstants.DBINWIDTH;
        }
    }    
    public double getCTermStart() { // y == c term
//System.out.println("in getCTermStart in DeChargedModifiedPeptideHit, " + mods.getDiffCTermMod());
        if(mods != null) {
            return ppl.getCTermStart() + mods.getDiffCTermMod();
        } else {     
            return ppl.getCTermStart();
        }
    }    
    public double getCTermDbinwidthStart() { // y == c term
        if(mods != null) {
            return (ppl.getCTermStart() + mods.getDiffCTermMod())* MassSpecConstants.DBINWIDTH;
        } else {

            return ppl.getCTermStart()* MassSpecConstants.DBINWIDTH;
        }
    }    

    // return sequence without leading and tailing residues
    public String getExactSequence() {
        if(exactSeq != null) {
            return exactSeq;
        }
        StringBuffer sb = new StringBuffer(getLength() + 15);
        int start = peptide.getStart();
        int end = peptide.getEnd();

        double nTermMod = mods == null? 0 : mods.getDiffNTermMod();
        double cTermMod = mods == null? 0 : mods.getDiffCTermMod();
        if(nTermMod != 0) {
            sb.append("(" + nTermMod + ")"); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != 0) {
            sb.append("(" + cTermMod + ")"); 
        }

        exactSeq = sb.toString(); 
        
        return exactSeq;
    }

    public String getExtendedSequence() {
        if(extendedSeq != null) {
            return extendedSeq;
        }
        StringBuffer sb = new StringBuffer(getLength() + 15);
        // 4 more bytes for the heading and tailing aa and '.'
        int len = peptide.getLength() + 4;
        int start = peptide.getStart();
        int end = peptide.getEnd();
        if (start == 0) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(start-1));
        }
        sb.append('.');

        TerminalModification nTermMod = mods.getNTerm();
        TerminalModification cTermMod = mods.getCTerm();
        if(nTermMod != null) {
            //sb.append("(" + nTermMod + ")"); 
            sb.append( nTermMod.getModInfo()); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != null) {
            //sb.append("(" + cTermMod + ")"); 
            sb.append(cTermMod.getModInfo()); 
        }
        sb.append('.');

        int lastIndex = peptide.getParent().getLength() - 1;
        if (end == lastIndex) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(end+1));
        }

        extendedSeq = sb.toString(); 
        
        return extendedSeq;
    }



    public void calcNumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        byte [] seq = getParent().getSequenceAsBytes();
        boolean [] boolMasses = ppl.getBoolMasses();
        MassCalculator mc = ppl.getMassCalculator();
        //fragMasses = mc.getFragMasses(); no big difference
        double yMass = getCTermStart(); // from super class
        double bMass = getNTermStart();  // from super class

        int start = peptide.getStart();
        int end = peptide.getEnd();
        int yIndex = start + end; // index for y ion, i is used for b ion

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi];
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);

            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int tempB = (int)(bMass * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            if(tempB < lastTrue && tempB > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempB]) numPeaksMatched++;
                //if(boolMasses[tempB+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(tempY < lastTrue && tempY > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempY]) numPeaksMatched++;

                //if(boolMasses[tempY+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());
    }

    public int [] getTheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
       
        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {
            int yi = yIndex -i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }
            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);
        }

        return theorMass;
    }
    public boolean equals(DeChargedModifiedPeptideHit m) {
        return super.equals(m) && mods.equals(m.mods);
    }
}
