/**
 * @file PeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.PeptideHit
 * @author Tao Xu
 * @date $Date: 2009/07/05 05:34:35 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.mspid.ProcessedPeakList;
import edu.scripps.pms.util.seq.Fasta;

import java.util.List;

public class PeptideHit implements Comparable<PeptideHit> {

    // description line of this PeptideHit sequence
    // the sequence string of this PeptideHit
    
    //private int chargeState;
    private int numPeaksMatched = 0;
    private int numTheroticPeaks = 0;
    private double probability;
    protected static double [] masses;
    protected ProcessedPeakList ppl;    
    protected Peptide peptide;
    public static final boolean isModified = false;

    public PeptideHit(Fasta parent, int start, int end) {
        peptide = new Peptide(parent, start, end);
    }  
    public PeptideHit(Peptide p) {
        peptide = p;
    }
    public void setProcessedPeakList(ProcessedPeakList ppl) {
        this.ppl = ppl;
        masses = ppl.getFragMasses();
//for(double f : masses) {
//    if(f != 0)
//    System.out.println(f + "\t");
//}
    }
    public boolean isModified() {
        return isModified;
    }
    public void calcNumPeaksMatched () {
        int mam = ppl.getSearchParams().getMultistageActivationMod();
  
        switch(ppl.getChargeState()) {
            case 3 : getCharge3NumPeaksMatched(mam); break;
            case 2 : getCharge2NumPeaksMatched(mam); break;  // both +1 and +2 considers only +1 fragments
            default : getChargeNNumPeaksMatched((ppl.getChargeState()+2)/2, mam); 
        }
    }

    public void libraryCalcNumPeaksMatched(int[] masses)
    {
        int mam = ppl.getSearchParams().getMultistageActivationMod();
        getNumPeaksMatched(mam,masses);
    }

    public double getTheorMass() {


        double mass = ppl.getMassCalculator().getPrecursorMass(peptide.getSequence()) +
                     ppl.getSearchParams().getStaticTerminalMods();
   
        return mass;
    }
    public int [] getTheorMasses() {
        // for multistage activation mode
        int mam = ppl.getSearchParams().getMultistageActivationMod();
        switch(ppl.getChargeState()) {
            case 3 : return getCharge3TheorMasses(mam); 
            case 2 : return getCharge2TheorMasses(mam);
            default : return getChargeNTheorMasses((ppl.getChargeState()+2)/2, mam); 
        }
        //return null;
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
    // multistage activation mod is only used by modified peptide hit
    public void getChargeNNumPeaksMatched(int z, int mam) {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        // y is used for both y and z, ie cterm ions
        // b is used for both b andc, ie nterm ions
        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function 
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;
        for(int i = start; i < end; i++) {
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);
//System.out.println((char)seq[i] + "\t" + bMass + "\t" + yMass);
            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/2f);
            // doubly charged fragment ions
            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;
 
            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
           // if(boolMasses[int2B++]) numPeaksMatched++;
           // if(boolMasses[int2Y++]) numPeaksMatched++;
            //numTheroticPeaks += 6;
            for(int j = 3; j <= z; j++) {
                int nb = (int)(0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/j);
                int ny = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/j);
                if(nb < lastTrue && nb > firstTrue) { 
                    numTheroticPeaks++;
                    if(boolMasses[nb]) numPeaksMatched++;
                }
                if(ny < lastTrue && ny > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[ny]) numPeaksMatched++;
                }
            }
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }

    public void getNumPeaksMatched(int mam, int [] masses) {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        // y is used for both y and z, ie cterm ions
        // b is used for both b andc, ie nterm ions
        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        for(int i = start; i < end; i++) {
            int mass = masses[i];
            if(mass < lastTrue && mass > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[mass]) numPeaksMatched++;
            }
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }

    // multistatge activation mod is only used in modified peptide hit function
    public void getCharge3NumPeaksMatched(int mam) { 
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        // y is used for both y and z, ie cterm ions
        // b is used for both b andc, ie nterm ions
        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function 
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        for(int i = start; i < end; i++) {
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);
//System.out.println((char)seq[i] + "\t" + bMass + "\t" + yMass);
            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/2f);
            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;
 
            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }
    public void getCharge2NumPeaksMatched(int mam) {
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
                //if(boolMasses[tempB+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(tempY < lastTrue && tempY > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[tempY]) numPeaksMatched++;

                //if(boolMasses[tempY+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
          // System.out.println("bmass: " + bMass + "\tymass: " + yMass); 
            //numTheroticPeaks += 4;
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());
    }
    // multistage activation is only used by modified peptide hit
    public int [] getChargeNTheorMasses(int z, int mam) {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;


        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
//if(bMass > 20)
//System.out.println("ystart: " + yMass + "\tbstart: " + bMass);
        //double tempb = 0; double tempy = 0; int indexb = 0; int indexy = 0;
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {
            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }
            bMass += masses[seq[i]];
            yMass += masses[seq[yIndex-i]];
            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);

            // for highly charged fragment ions
            processHighlyChargedIons(theorMass, bMass, yMass, z);
        }
        return theorMass;
    }
    // multistage activation mod is only used in modified peptide hits
    public int [] getCharge3TheorMasses(int mam) {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();


//if(bMass > 20)
//System.out.println("ystart: " + yMass + "\tbstart: " + bMass);
        //double tempb = 0; double tempy = 0; int indexb = 0; int indexy = 0;
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            bMass += masses[seq[i]];
            yMass += masses[seq[yIndex-i]];
            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);
        }
        return theorMass;
    }

    public int [] getCharge2TheorMasses(int mam) {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;


        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
        
        //double tempb = 0;
        //double tempy = 0;
//if(bMass > 20)
//System.out.println("ystart: " + yMass + "\tbstart: " + bMass);

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }
            bMass += masses[seq[i]];
            yMass += masses[seq[yIndex-i]];
            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);
        }
        return theorMass;
    }

    // n-term ion: b ion for cid and c ion for etd
    protected void processSinglyChargedBIon(int [] theorMass, double mass) {
//System.out.println("in singlyChargedBIon, mass: " + mass);
        try {
            
            assignTheorMass(theorMass, mass, 1);
            //int intMass = (int)mass;
            //theorMass[intMass] = 50;
            //int index = intMass + 1;
            //if(theorMass[index] < 25) theorMass[index] = 25;

            //index++; // -= 2; 
            //if(theorMass[index] < 25) theorMass[index] = 25;

            if(!ppl.isEtdSpectrum()) {
                //index = intMass - 18; // H2O loss
                //if(theorMass[index] < 10) theorMass[index] = 10;
                assignValue(theorMass, mass-18, 10);
                //index = intMass - 28; // CO loss
                //if(theorMass[index] < 10) theorMass[index] = 10;
                assignValue(theorMass, mass-28, 10);
                //index = intMass - 17; // NH3 loss
                //if (theorMass[index] < 10) theorMass[index] = 10;
                assignValue(theorMass, mass-17, 10);
            }
        } catch(Exception e) {}// igore exception caused by weird aa residue
       
    }
    // c-term ion: y ion for cid and z ion for etd
    protected void processSinglyChargedYIon(int [] theorMass, double mass) {
//System.out.println("in singlyChargedYIon, mass: " + mass);
        try { 
            assignTheorMass(theorMass, mass, 1);
            //int intMass = (int)mass;
            //System.out.println("intmass: " + intMass);
            //theorMass[intMass] = 50;
            //int index = intMass + 1;
            //if(theorMass[index] < 25)  theorMass[index] = 25;
 
            //index++; // -= 2; 
            //if(theorMass[index] < 25) theorMass[index] = 25;

            if(!ppl.isEtdSpectrum()) {
                //index = intMass - 17; // loss NH3
                //if(theorMass[index] < 10) theorMass[index] = 10;
                assignValue(theorMass, mass-17, 10);
            }
        
        } catch(Exception e) {}
    }
    // c-term ion: y ion for cid and z ion for etd
    protected void processHighlyChargedIons(int theorMass[], double bMass, double yMass, int z) {
//System.out.println("in doublyChargedYIon, mass: " + yMass);
        try {
            for(int i = 3; i <= z; i++) {
                assignTheorMass(theorMass, bMass + (i-1)*MassSpecConstants.MASSPROTONDB, i);
                assignTheorMass(theorMass, yMass + (i-1)*MassSpecConstants.MASSPROTONDB, i);
            }
        } catch(Exception e) {}// igore exception caused by weird aa residue
    }
    protected void assignTheorMass(int theorMass[], double mass, int z) {
        int averagineIndex = (int)mass/500;

//System.out.println("mass: " + mass + "\tz: " + z + "\tindex: " + averagineIndex);
        switch(averagineIndex) {
            case 0 : 
                     assignValue(theorMass, mass/z, 50);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 14);
                     break;
            case 1 : 
                     assignValue(theorMass, mass/z, 50);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 28);
                     break;
            case 2 : 
                     assignValue(theorMass, mass/z, 50);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 43);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 22);
                     break;
            case 3 : 
                     assignValue(theorMass, mass/z, 45);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 30);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 15);
                     break;
            case 4 : 
                     assignValue(theorMass, mass/z, 36);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 39);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 22);
                     break;
            case 5 : 
                     assignValue(theorMass, mass/z, 30);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 46);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 30);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 15);
                     break;
            case 6 : 
                     assignValue(theorMass, mass/z, 24);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 47);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 37);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 22);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 11);
                     break;
            case 7 : 
                     assignValue(theorMass, mass/z, 19);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 42);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 43);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 29);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 15);
                     break;
            case 8 : 
                     assignValue(theorMass, mass/z, 15);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 38);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 47);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 35);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 21);
                     assignValue(theorMass, (mass + 6*MassSpecConstants.MASSDIFFC12C13)/z, 11);
                     break;
            case 9 : 
                     assignValue(theorMass, mass/z, 12);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 33);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 49);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 40);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 27);
                     assignValue(theorMass, (mass + 6*MassSpecConstants.MASSDIFFC12C13)/z, 15);
                     break;
            case 10 : 
                     assignValue(theorMass, mass/z, 9);
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 28);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 45);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 44);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 31);
                     assignValue(theorMass, (mass + 6*MassSpecConstants.MASSDIFFC12C13)/z, 19);
                     assignValue(theorMass, (mass + 7*MassSpecConstants.MASSDIFFC12C13)/z, 11);
                     break;
            case 11 : 
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 24);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 42);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 47);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 36);
                     assignValue(theorMass, (mass + 6*MassSpecConstants.MASSDIFFC12C13)/z, 24);
                     assignValue(theorMass, (mass + 7*MassSpecConstants.MASSDIFFC12C13)/z, 14);
                     break;
            case 12 : 
                     assignValue(theorMass, (mass + MassSpecConstants.MASSDIFFC12C13)/z, 20);
                     assignValue(theorMass, (mass + 2*MassSpecConstants.MASSDIFFC12C13)/z, 38);
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 42);
                     assignValue(theorMass, (mass + 6*MassSpecConstants.MASSDIFFC12C13)/z, 29);
                     assignValue(theorMass, (mass + 7*MassSpecConstants.MASSDIFFC12C13)/z, 18);
                     break;
            default : 
                     assignValue(theorMass, (mass + 3*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 4*MassSpecConstants.MASSDIFFC12C13)/z, 50);
                     assignValue(theorMass, (mass + 5*MassSpecConstants.MASSDIFFC12C13)/z, 42);
                     break;
        }
    } 
    protected void assignValue(int theorMass[], double mass, int minValue) {
        //int index = (int)(mass + 0.5);
        if(mass < ppl.getPeakList().getMinM2z() || mass > ppl.getPeakList().getMaxM2z()) {
            return;
        }
        int index = (int)(mass);
        if(theorMass[index] < minValue) {
            theorMass[index] = minValue;
        }
    }
    // c-term ion: y ion for cid and z ion for etd
    protected void processDoublyChargedYIon(int theorMass[], double yMass) {
//System.out.println("in doublyChargedYIon, mass: " + yMass);
        try {
            assignTheorMass(theorMass, yMass+MassSpecConstants.MASSPROTONDB+0.5, 2);
            double tempy = (yMass+MassSpecConstants.MASSPROTONDB)/2.f;
            //int indexY = (int)(tempy + 0.25);
            //theorMass[indexY] = 50;
            //int index = indexY + 1;
            //if(theorMass[index] < 25) theorMass[index] = 25;

            //index++; // -= 2; // may not be right
            //if(theorMass[index] < 25) theorMass[index] = 25;

            if(!ppl.isEtdSpectrum()) {
                //index = (int)(tempy - ProcessedPeakList.PLUS3NH3IONADJUSTMENT); 
                assignValue(theorMass, tempy - ProcessedPeakList.PLUS3NH3IONADJUSTMENT, 10);
                //if(theorMass[index] < 10) theorMass[index] = 10; 
            }
        } catch(Exception e) {}// igore exception caused by weird aa residue
    }

    // n-term ion: b ion for cid and c ion for etd
    protected void processDoublyChargedBIon(int theorMass[], double bMass) {
//System.out.println("in doublyChargedBIon, mass: " + bMass);
        try {
            assignTheorMass(theorMass, bMass+MassSpecConstants.MASSPROTONDB+0.5, 2);
            double tempb = (bMass + MassSpecConstants.MASSPROTONDB)/2.f;
            //int indexB = (int)(tempb + 0.25);
            //theorMass[indexB] = 50;
            //int index = indexB + 1;
            //if(theorMass[index] < 25) theorMass[index] = 25;

            //index++; // -= 2; // may not be right
            //if(theorMass[index] < 25) theorMass[index] = 25;
            if(!ppl.isEtdSpectrum()) {
                //index = indexB - 9; // for loss H2O
                assignValue(theorMass, tempb-9, 10);
                //if(theorMass[index] < 10) theorMass[index] = 10;
                //index = indexB - 14; // for loss CO
                //if(theorMass[index] < 10) theorMass[index] = 10;
                assignValue(theorMass, tempb-14, 10);
              
 
                // may not be accurate 
                //index = (int)(tempb - ProcessedPeakList.PLUS3NH3IONADJUSTMENT); // for loss NH3
                //if(theorMass[index] < 10) theorMass[index] = 10;
                assignValue(theorMass, tempb - ProcessedPeakList.PLUS3NH3IONADJUSTMENT, 10);
            }
        
        } catch(Exception e) {}// igore exception caused by weird aa residue
    }
    public boolean equals(PeptideHit p) {
        
        if(p != null) {
            return getSequence().equals(p.getSequence()); 
        }
        return false;
    }
    public int compareTo(PeptideHit p) {

//        return p.numPeaksMatched - this.numPeaksMatched; 
        if(p.probability < this.probability) {
            return 1;
        } else if(p.probability > this.probability) {
            return -1;
        } else { 
            if(p.isModified() == this.isModified()) {
                return 0; 
            } else if(p.isModified()) {
                return -1;
            } else {
                return 1;
            }
        }
       
    }
    public void setNumPeaksMatched(int n, int t, double ptrue) {
        numPeaksMatched = n;
        numTheroticPeaks = t;
        //probability = DistributionCalculator.getBinomialSum((int)(ptrue*100+0.5), t, n);
       probability = DistributionCalculator.getBinomialSum(ptrue, t, n);
    }
    public Fasta getParent() {
        return peptide.getParent();
    }
    public double getProbability() {
        return probability;
    }
    public String getExtendedSequence() {
        return peptide.getExtendedSequence();
    }
    public String getExtraExtendedSequence() {
        return peptide.getExtraExtendedSequence();
    }
    public String getSequence() {
        return peptide.getSequence(); 
    }
    public String getExactSequence() {
        return peptide.getSequence(); 
    }
    public int getStart() {
        return peptide.getStart();
    }
    public String getDefline() {
        return peptide.getDefline();
    }
    public int getEnd() {
        return peptide.getEnd();
    }

    public int getNumPeaksMatched() {
        return numPeaksMatched;
    }
    public int getLength() {
        return peptide.getLength();
    }
    public String getSequestLikeAccession() {
        return peptide.getSequestLikeAccession();
    }
    public List<String> getDefList() {

        return peptide.getDefList();
    }
    public String getAccession() {
        return peptide.getAccession();
    }
    public int getNumPeaks() {
        //return 2*2*(chargeState-1)*(end - start); // 2*(lengh - 1)
        return numTheroticPeaks;
    }
    public DiffMod [] getDiffMods()
    {
        return null;
    }

    public Modifications getModifications()
    {
        return null;
    }
}
