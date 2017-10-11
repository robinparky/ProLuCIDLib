/**
 * @file DeChargedPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.DeChargedPeptideHit
 * @author Tao Xu
 * @date $Date: 2009/02/02 05:11:43 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;

public class DeChargedPeptideHit extends PeptideHit {

    // description line of this DeChargedPeptideHit sequence
    // the sequence string of this DeChargedPeptideHit
    
    //private int chargeState;

    public DeChargedPeptideHit(Fasta parent, int start, int end) {
        super(parent, start, end); 
    }  
    public DeChargedPeptideHit(Peptide p) {
        super(p);
    }
    public void calcNumPeaksMatched () {

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
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int tempB = (int)(bMass * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            if(tempB < lastTrue && tempB > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[tempB]) numPeaksMatched++;
            }
            if(tempY < lastTrue && tempY > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[tempY]) numPeaksMatched++;

            }
            
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());
    }
    public double getTheorMass() {
        double mass = ppl.getMassCalculator().getPrecursorMass(peptide.getSequence()) +
                     ppl.getSearchParams().getStaticTerminalMods();
   
        return mass;
    }
    public int [] getTheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;
        //double tempb = 0;
        //double tempy = 0;
//if(bMass > 20)
//System.out.println("ystart: " + yMass + "\tbstart: " + bMass);

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            bMass += masses[seq[i]];
            yMass += masses[seq[yIndex-i]];
            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);
        }
        return theorMass;
    }
    
    public double getNTermDbinwidthStart() {
        return ppl.getNTermDbinwidthStart();
    }
    public double getCTermDbinwidthStart() {
        return ppl.getCTermDbinwidthStart();
    }
    public double getNTermStart() {
        return ppl.getNTermStart();
    }
    public double getCTermStart() {
        return ppl.getCTermStart();
    }


    // n-term ion: b ion for cid and c ion for etd
    protected void processSinglyChargedBIon(int [] theorMass, double mass) {
//System.out.println("in singlyChargedBIon, mass: " + mass);
        try {
            assignTheorMass(theorMass, mass);
            int intMass = (int)mass;
            if(!ppl.isEtdSpectrum()) {
                int index = intMass - 18; // H2O loss
                if(theorMass[index] < 10) theorMass[index] = 10;
                index = intMass - 28; // CO loss
                if(theorMass[index] < 10) theorMass[index] = 10;
                index = intMass - 17; // NH3 loss
                if (theorMass[index] < 10) theorMass[index] = 10;
            }
        } catch(Exception e) {}// igore exception caused by weird aa residue
       
    }
    // c-term ion: y ion for cid and z ion for etd
    protected void processSinglyChargedYIon(int [] theorMass, double mass) {
//System.out.println("in singlyChargedYIon, mass: " + mass);
        try { 
            assignTheorMass(theorMass, mass, 1);
            int intMass = (int)mass;
            if(!ppl.isEtdSpectrum()) {
                int index = intMass - 17; // loss NH3
                if(theorMass[index] < 10) theorMass[index] = 10;
            }
        
        } catch(Exception e) {}
    }
    protected void assignTheorMass(int theorMass[], double mass) {

         assignValue(theorMass, mass, 50);
    } 
}
