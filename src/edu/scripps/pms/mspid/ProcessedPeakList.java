/**
 * @file ProcessedPeakList.java
 * This is the source file for edu.scripps.pms.mspid.ProcessedPeakList
 * @author Tao Xu
 * @date $Date: 2009/09/01 05:24:59 $
 */
package edu.scripps.pms.mspid;

import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.spectrum.*;

import java.io.*;
import java.util.*;

public class ProcessedPeakList {
    //int numGreater = 0;
    public static final double PLUS3NH3IONADJUSTMENT =  8.5* MassSpecConstants.DBINWIDTH - 0.25;
    public static final int DEFAULTPREPROCESS = 0;
    public static final int XCORRPREPROCESS = 1;
    public static final int TOPDOWNPREPROCESS = 2;
    public static final int XCORRWINDOWSIZE = 1100; 
    //public static final int MAXNUMPEAKS = 80000;
    public static final int PRECISIONFACTOR = 1000; // for sp process
    public static final int ACCURACYFACTOR = 1;  // for XCorr
    
    protected PeakList peakList;
    protected ArrayList<Peak> peaks;
    protected SearchParams params;
    //protected double [] intensityVal = new double[MAXNUMPEAKS];
    protected double [] intensityVal; // = new double[MAXNUMPEAKS];
    protected double [] intensityVal2;
    
    // as exp_intensity in pep_prob
    //protected boolean [] boolMasses = new boolean[MAXNUMPEAKS]; // boolean representation of masses
    protected boolean [] boolMasses; // = new boolean[MAXNUMPEAKS]; // boolean representation of masses
//    protected int numPeaks; 
    protected double cutoff = 0.06;
    //protected int lowestM2z; 
    //protected int highestM2z;
    protected double [] massCorr; // stores the transformed intensity
    protected double [] massSum; // stores the sum of +- 75 intensities to speed up the correlation
    protected double prcMass; // this is real precursor M+H from z line
    protected int chargeState;
    protected boolean isDeCharged = false;
    //double tolerance;
    double fragTolerance;
    //double highLimit;
    //double lowLimit; 
    protected double entropy = 0;
    protected double totalIntensity = 0;
    protected double totalTransformedIntensity = 0;
    protected double maxTransformedInten = 0;
    protected Zline zline;
    //protected double [] tempIntensity; 
    protected double globMax; // the cutoff for transformed intensity
    protected int maxShift;
    double maxShiftProduct;
    protected int minimumPeptideLength;
    protected int numTrues = 0;
    protected int firstTrue = 0;
    protected int lastTrue = 0;
    protected double pTrue;  // propotion of trues in boolMasses
    protected int numFragBins; // number of bins that could be true
    protected int finalNumPeaks; // final number of peaks

    protected MassCalculator mc;
    protected double [] masses; // DBINWIDTH masses
//    protected int lowM2z;
//    protected int highM2z;

    //protected Modifications mods;

    protected double cTermStart; // c term for y and z ions
    protected double nTermStart; // n term for b and c ions

    protected double cTermDbinwidthStart;
    protected double nTermDbinwidthStart;

    protected boolean isEtd;
//    protected List<SearchResult> searchResults = new ArrayList();

//    public abstract int calcNumPeaksMatched (PeptideHit p);
    //public abstract int [] calcNumPeaksMatched(Fasta f, int start, int end);

//    protected abstract int [] getTheorMasses(ScoredPeptideHit s);

    public boolean isEtdSpectrum() {
        return isEtd; 
    }
    public double [] getFragMasses() {
        return masses;
    } 
    public double getNTermStart() {
        return nTermStart;
    }
    public double getCTermStart() {
        return cTermStart;
    }
    public double getCTermDbinwidthStart() {
        return cTermDbinwidthStart;
    }
    public double getNTermDbinwidthStart() {
        return nTermDbinwidthStart;
    }
    public int getChargeState() {
        return chargeState;
    }
    public Iterator<Modifications> getModifications() {
        return params.getAllModifications();
    }


    public ProcessedPeakList(PeakList peaklist, Zline z, SearchParams sp, MassCalculator mc, boolean preprocess) {
        init(peaklist,z,sp,mc,preprocess);
    }

    private void init(PeakList peaklist, Zline z, SearchParams sp, MassCalculator mc, boolean preprocess)
    {
        this.peakList = peaklist;
        this.mc = mc;
        zline = z;
        prcMass = z.getM2z();// precursor mass
        chargeState = z.getChargeState() > 2? z.getChargeState() : 2;
        isEtd = peakList.getActivationType() != null && peakList.isEtdSpectrum() && sp.isEtdSearch();
        //System.out.println("Scan#: " + peaks.getLoscan() + "\t+" + z.getChargeState() + "\tNumPeaks: " + peaks.numPeaks()+ "\tm2z: " + z.getM2z() + "ptrue: " + pTrue);
        params = sp;
        isDeCharged = params.isDeCharged();
        masses = mc.getFragMasses(MassSpecConstants.DBINWIDTH);
        //       lowM2z = (int)((peaks.getMinM2z()-10)*MassSpecConstants.DBINWIDTH + 0.5f);
        //       highM2z = (int)((peaks.getMaxM2z()+10)*MassSpecConstants.DBINWIDTH + 0.5f);

        double yStart = MassSpecConstants.MASSH3O + params.getStaticCTermMod();
        double bStart = MassSpecConstants.MASSPROTON + params.getStaticNTermMod();

        cTermStart = isEtd? yStart-16.0187 : yStart;
        nTermStart = isEtd? bStart+17.02655 : bStart;

        //System.out.println("in ProcessedPeakList, isEtd: " + isEtd + "\tcTermStart: " + cTermStart + "\tnTermStart: " + nTermStart);
        cTermDbinwidthStart = cTermStart* MassSpecConstants.DBINWIDTH;
        nTermDbinwidthStart = nTermStart* MassSpecConstants.DBINWIDTH;
//        System.out.println("in ProcessedPeakList, isEtd: " + isEtd + "\tcTermDbinStart: " + cTermDbinwidthStart + "\tnTermDbinwidthStart: " + nTermDbinwidthStart);

    //    massCorr = new double[(int)(prcMass* MassSpecConstants.DBINWIDTH + params.getMaxMassShift() + 20.5f)*ACCURACYFACTOR];
   //     massSum = new double[massCorr.length];

        this.peaks = new ArrayList<Peak>(peakList.numPeaks());
        //boolMasses = new boolean[(int) (prcMass+params.getMaxMassShift()+20.5)*PRECISIONFACTOR]; // = new boolean[MAXNUMPEAKS]; // boolean representation of masses
//        intensityVal = new double [(int)(prcMass+params.getMaxMassShift()+20.5)*PRECISIONFACTOR + 100]; // = new double[MAXNUMPEAKS];

        fragTolerance = params.getFragmentTolerance()/1000.0f; // change from ppm to .4 etc.
        minimumPeptideLength = params.getMinimumPeptideLength() - 2; // to avoid +1 and >=

        //highLimit = prcMass + params.getHighPrecursorTolerance();
        //lowLimit = prcMass - params.getLowPrecursorTolerance();

        maxShift = (int)(PRECISIONFACTOR*fragTolerance);
//System.out.println("maxshift: " + maxShift + "\tfragTolerance : " + fragTolerance + "\tboolmass length: " + boolMasses.length);
        maxShiftProduct = maxShift*maxShift;
        maxTransformedInten = 0;
        for(Iterator<Peak> peakIt = peakList.getPeaks(); peakIt.hasNext();) {
            Peak p = peakIt.next();
            if(p.getM2z() < prcMass) {
                peaks.add(p);
            }
            double inten = p.getIntensity(); // intensity
            //tempIntensity[counter] = Math.sqrt(inten);
            double sqrtIntens = Math.sqrt(inten);
            totalIntensity += inten;  // calculate the totalIntensity
            totalTransformedIntensity += sqrtIntens;
            if (sqrtIntens > maxTransformedInten) {
                maxTransformedInten = sqrtIntens;
//System.out.println("Max intens: " + inten + "\tmaxTransformedIntens: " + maxTransformedInten);
            }
        }
        globMax = maxTransformedInten * cutoff;
        globMax = globMax > 2? globMax : 2;  // this does not work well with orbtrap data
        if(preprocess){
            preprocess(sp.getPreprocess());
    //        buildExpMass();
        }
    }

    public double getPrecursorMass() {
        return prcMass;
    }
    public boolean isDeCharged() {
        return isDeCharged;
    }
    protected void processSinglyChargedBIon(int [] theorMass, double mass) {
        try {
            int intMass = (int)mass;
            theorMass[intMass] = 50;
            int index = intMass + 1;
            if(theorMass[index] < 25) theorMass[index] = 25;

            index++; // -= 2; 
            if(theorMass[index] < 25) theorMass[index] = 25;

            index = intMass - 18; // H2O loss
            if(theorMass[index] < 10) theorMass[index] = 10;
            index = intMass - 28; // CO loss
            if(theorMass[index] < 10) theorMass[index] = 10;
            index = intMass - 17; // NH3 loss
            if (theorMass[index] < 10) theorMass[index] = 10;
        } catch(Exception e) {}// igore exception caused by weird aa residue
       
    }
    protected void processSinglyChargedYIon(int [] theorMass, double mass) {
        try { 
            int intMass = (int)mass;
            //System.out.println("intmass: " + intMass);
            theorMass[intMass] = 50;
            int index = intMass + 1;
            if(theorMass[index] < 25)  theorMass[index] = 25;
 
            index++; // -= 2; 
            if(theorMass[index] < 25) theorMass[index] = 25;

            index = intMass - 17; // loss NH3
            if(theorMass[index] < 10) theorMass[index] = 10;
        
        } catch(Exception e) {}
    }
    protected void processDoublyChargedYIon(int theorMass[], double yMass) {
        try {
            double tempy = (yMass+ MassSpecConstants.MASSPROTONDB)/2.f;
            int indexY = (int)(tempy + 0.25);
            theorMass[indexY] = 50;
            int index = indexY + 1;
            if(theorMass[index] < 25) theorMass[index] = 25;

            index++; // -= 2;
            if(theorMass[index] < 25) theorMass[index] = 25;

            index = (int)(tempy - PLUS3NH3IONADJUSTMENT); 
            if(theorMass[index] < 10) theorMass[index] = 10; 
        } catch(Exception e) {}// igore exception caused by weird aa residue?
    }

    protected void processDoublyChargedBIon(int theorMass[], double bMass) {
        try {
            double tempb = (bMass + MassSpecConstants.MASSPROTONDB)/2.f;
            int indexB = (int)(tempb + 0.25);
            theorMass[indexB] = 50;
            int index = indexB + 1;
            if(theorMass[index] < 25) theorMass[index] = 25;

            index++; // -= 2;
            if(theorMass[index] < 25) theorMass[index] = 25;

            index = indexB - 9; // for loss H2O
            if(theorMass[index] < 10) theorMass[index] = 10;
            index = indexB - 14; // for loss CO
            if(theorMass[index] < 10) theorMass[index] = 10;
            index = (int)(tempb - PLUS3NH3IONADJUSTMENT); // for loss NH3
            if(theorMass[index] < 10) theorMass[index] = 10;
        
        } catch(Exception e) {}// igore exception caused by weird aa residue
    }
    public double autoCorrelation(int [] values) {
        double sumProduct = 0.f;
        double moreSumProduct= 0.f;
        int numElements = values.length;
        for(int i=0; i < numElements; i++) {
            sumProduct += values[i]*values[i];
        }
        for(int i=0; i < numElements; i++) {
            if(values[i] != 0.f)
                for(int j=1; j <=75; j++)
                    if(i+j < numElements)
                        if(values[i+j] != 0.f)
                            moreSumProduct += values[i]*values[i+j];
        }

        moreSumProduct *= 2;
        return sumProduct - (moreSumProduct)/151.f;
    }
    
    //public abstract void correlation(ScoredPeptideHit s);
    public void correlation(ScoredPeptideHit s) {
        int theorMass[] = s.getTheorMasses();
       // System.out.println(s.getSequence() + " cs "+this.getChargeState());

      /*  if(s.isModified())
        {
            DiffMod [] diffs = s.getDiffMods();
            Modifications mods = s.getModifications();
            double cterm = mods.getDiffCTermMod();
            double nterm = mods.getDiffNTermMod();
            System.out.println("CTerm is "+cterm);
            System.out.println("Nterm is "+nterm);
            if(diffs!=null)
            {
                System.out.println("DbinwidthMassShift");

                for(DiffMod diff: diffs)
                {
                    if(diff==null ) System.out.print(" 0 ");
                    else System.out.print(" " +diff.getDbinwidthMassShift() +" ");

                }
                System.out.println();
                System.out.println("DbinwidthNeutralLoss");
                for(DiffMod diff: diffs)
                {

                    if(diff==null ) System.out.print(" 0 ");
                    else System.out.print(" "+diff.getDbinwidthNeutralLoss()+" ");
                }
                System.out.println();
            }
        }*/
       /* TheoreticalPeaks tp = new TheoreticalPeaks(getSearchParams(),getMassCalculator());
        int theorMass2[] = tp.calcTheoreticalPeaks(this,s.getSequence().getBytes(),s.getPeptideHits().get(0).isModified());
        System.out.println(Arrays.toString(theorMass));
        System.out.println(Arrays.toString(theorMass2));
        for(int i=0; i<theorMass.length;i++)
        {
            if(theorMass[i]!=theorMass2[i])
            {
                System.out.println("not same");
                break;
            }
        }*/
        ArrayList<Peak> theorPeaks = new ArrayList<Peak>(200);

        int maxIndex = massCorr.length - 75;
        for(int i = 0; i < theorMass.length; i++) {
            if(theorMass[i] > 0 && i > 75 && i < maxIndex) {
                theorPeaks.add(new Peak(i, theorMass[i]));
            }
        }



/*
System.out.println("massCorr");
for(int i = 0; i < massCorr.length; i++) {
    if(massCorr[i] > 0) {
        System.out.println(i + "\t" + massCorr[i]);
    }
}
System.out.println("thereomasses");
for(Peak p : theorPeaks) { 
    int mass = (int)p.getM2z()*ACCURACYFACTOR;
    double intensity = p.getIntensity();

    System.out.println(mass + "\t" + intensity);
}
*/
        double sumProduct = 0;
        double moreSumProduct = 0; 
        double newMoreSumProduct = 0; 
        for(Peak p : theorPeaks) {
            int mass = (int)p.getM2z()*ACCURACYFACTOR;
            double intensity = p.getIntensity();
            sumProduct += massCorr[mass]*intensity;
            newMoreSumProduct += massSum[mass]*intensity;
//            for(int i = 1; i <= 75; i++) {
//                moreSumProduct += intensity*massCorr[mass+i]; 
//                moreSumProduct += intensity*massCorr[mass-i]; 
//            }
           
        }

        //double xcorr = (0.993377483f*sumProduct - moreSumProduct/151)/10000;
        double xcorr = (0.993377483f*sumProduct - newMoreSumProduct/151)/10000;
        xcorr = xcorr < 0.00001? 0.00001f : xcorr;

 
        s.setXCorr(xcorr);
        //s.setScore(xcorr);
    }

    public ScoredPeptideHit prelimScoreCorrelation(ProcessedPeakList peakList)
    {
        int numTheroticPeaks =0;
        int numPeaksMatched =0;
        int i=0;
        boolean [] boolMass = peakList.getBoolMasses();
        boolean [] mybools = this.getBoolMasses();
        for(int j=firstTrue; j<lastTrue; j++)
        {
            if(boolMass[j] )
            {
                numTheroticPeaks++;
              if(mybools[j]) numPeaksMatched++;
            }
        }
        peakList.dumpBoolMass();
        //this.dumpBoolMass();
        double probability = DistributionCalculator.getBinomialSum(this.getPTrue(), numTheroticPeaks, numPeaksMatched);
        ScoredPeptideHit sph = new ScoredPeptideHit(probability);
        return sph;
    }



    public ScoredPeptideHit prelimScoreCorrelation(PeakList peakList)
    {
        int numTheroticPeaks =0;
        List<Peak> sortedPeak = peakList.getSortedPeaks(PeakList.SORTBYINTENSITY);
        int numPeaksMatched =0;
        int i=0;
        Zline z= peakList.getZlines().next();
        for (Peak peak: sortedPeak )
        {
            double m2z = peak.getM2z();
            if(peak.getM2z()<z.getM2z()) {
                int intM2z = (int) (PRECISIONFACTOR*m2z + 0.5f);

                //if (sqrtIntens > globMax) {
                int indexShift = maxShift;
                while (indexShift >= 0) {
                    //double maxShiftProduct = maxShift*maxShift;
                    // transformedInt as ftemp in PREPROCESS in pep_prob
                    int left = intM2z-indexShift;
                    int right = intM2z+indexShift;
                    if(intM2z>firstTrue && intM2z<lastTrue)
                    {
                        numTheroticPeaks++;
                        //numTheroticPeaks++;
                        if(boolMasses[left])
                        {
                            numPeaksMatched++;
                        }
                        if(boolMasses[right])
                        {
                            numPeaksMatched++;
                            if(boolMasses[left])numTheroticPeaks++;
                        }
                    }

                    indexShift--;
                }
            }


            i++;
           // if(i>=sortedPeak.size()/2) break;
        }
        double probability = DistributionCalculator.getBinomialSum(this.getPTrue(), numTheroticPeaks, numPeaksMatched);
        ScoredPeptideHit sph = new ScoredPeptideHit(probability);
        return sph;
    }


    public MassCalculator getMassCalculator() {
        return mc;
    }
    public SearchParams getSearchParams() {
        return params;
    }
    public Zline getZline() {
        return zline;
    }
    public PeakList getPeakList() {
        return peakList;
    }

    public boolean [] getProcessedMasses() {
        return boolMasses;
    }
    public double getEntropy() {
        if (entropy == 0) {
            entropy = calcEntropy();
        }
        return entropy;
    }
    // calculate entropy
    protected double calcEntropy() {
        double sum = 0;
        for(Peak p : peaks) {
            double ratio = p.getIntensity()/totalIntensity;
            sum += -ratio* Math.log(ratio);
        }
        return sum/(double) Math.log(2.0f);

    }
    public double [] getMassCorr() {
        return massCorr;
    }
    protected boolean isInPrecursorRange(ArrayList<Range> ranges, double mz) {
        for(Iterator<Range> it = ranges.iterator(); it.hasNext();) {
            Range r = it.next();
            if(r.isInRange(mz)) {
                return true;
            }
        }
        return false;
    }
    public void buildExpMass() {
        //double prcMass = peaks.getPrecursorMass(); // precursor massa
        double maxMass = prcMass*ACCURACYFACTOR; // zline.getM2z(); //peaks.getPrecursorMass();//getMaxM2z() + 50; // ftemp in BUILD_EXP in pepprobe
        int highestIndex = 0;
        double maxIntensity = 0f;
        ArrayList<Range> precursorRanges = peakList.getPrecursorRanges(zline);
//for(Range r : precursorRanges) {
//System.out.println("Low: " + r.getLowBound() + "\thigh " + r.getHighBound());
//}
        for(Peak p : peaks) {
            if(isInPrecursorRange(precursorRanges, p.getM2z())) {
//System.out.println(p.getM2z() + "\t intensity: " + p.getIntensity());
                continue;
            }
            double m2zF = p.getM2z()*ACCURACYFACTOR;
            if(m2zF < prcMass*ACCURACYFACTOR) {
                int m2z = (int)(m2zF* MassSpecConstants.DBINWIDTH + 0.5f);
                double intensity = (double) Math.sqrt((int)p.getIntensity());
                if(massCorr[m2z] < intensity) {
                    massCorr[m2z] = intensity;
                    if(intensity > maxIntensity) {
                        maxIntensity = intensity;
                    }
                }
            }
            
        }
        double overAll = 0f;
        double tempIntens[] = new double[massCorr.length]; // mass_temp
        
        for(int i = 0; i < massCorr.length; i++) {
            if(massCorr[i] != 0) {
                tempIntens[i] = massCorr[i]*100.f/maxIntensity; // adjust all intensity to [0, 100]
                if(tempIntens[i] > overAll) {
                    overAll = tempIntens[i]; // junck code? overAll == 100??
                }
                massCorr[i] = 0;
            }
        }
        
        int win = 10;
        //int winSize = (int)(peaks.getMaxM2z()*MassSpecConstants.DBINWIDTH + 0.5)/win; 
        int winSize = (int)(prcMass* MassSpecConstants.DBINWIDTH*ACCURACYFACTOR + 0.5)/win;
        overAll *= 0.05f; // make overAll == 5

        for(int i = 0; i < win; i++) {
            double max = 0;
            for(int j = 0; j < winSize; j++) {
                int index = i*winSize + j;
                if(tempIntens[index] > max) {
                    max = tempIntens[index];
                }
            }
            for(int j = 0; j < winSize; j++) {
                int index = i*winSize + j;
                if(tempIntens[index] > overAll) {
                    //peaks with tempIntens < 5 will be ignored, i.e, peaks 
                    // with original intensity <  25 will be ignored
                    massCorr[index] = tempIntens[index]*50/max; // make all peaks to [2.5, 50]
                }
            }
        }
       /* 
        int numNone0 = 0;
        double min = 1000;
        for(int i = 0; i< massCorr.length; i++) {
          
            if(massCorr[i] > 0) {
                System.out.println(i + "\t" + massCorr[i]);
                if(massCorr[i] < min) {
                    min = massCorr[i];
                }
                numNone0++;
            }
        }   
        System.out.println("Num None 0 massCorr peaks: " + numNone0 + "\tmin: " + min);
       */

        int numBins = massCorr.length - 1;
        for(int i = 0; i < numBins; i++) {
            int left = i - 75;
            int right = i + 75;
            left = left < 0 ? 0 : left;
            right = right > numBins ? numBins : right;
            double intensitysum = 0;
            for(int j = left; j < i; j++) {
                intensitysum += massCorr[j];
            }
            for(int j = right; j > i; j--) {
                intensitysum += massCorr[j];
            }
            massSum[i] = intensitysum;
        }

    }
    public void preprocess(int mode) {
        TimeUtils timer = new TimeUtils();
        //timer.startTiming();
        switch (mode) {
            case DEFAULTPREPROCESS: defaultPreprocess();
            break;
            case XCORRPREPROCESS: xcorrPreprocess();
               // writeLogFile();
            break;
            case TOPDOWNPREPROCESS: topdownPreprocess();
            break;
        }
        //System.out.print("Scan#: " + peakList.getLoscan() + "\t+" + zline.getChargeState() + "\tNumPeaks: " + peakList.numPeaks()+ "\tm2z: " + zline.getM2z() + "\tptrue: " + pTrue);
    }

   
    protected void shift(int indexShift, int intM2z, double transformedIntensity) {
        
        int leftIndex = intM2z - indexShift;
        int rightIndex = intM2z + indexShift;
        boolMasses[leftIndex] = true;
        boolMasses[rightIndex] = true;
/*        if (intensityVal[leftIndex] <= transformedIntensity) {
            intensityVal[leftIndex] = transformedIntensity;
        } 
           
        if (intensityVal[rightIndex] <= transformedIntensity) {
            intensityVal[rightIndex] = transformedIntensity;
        }*/
    }
    protected void topdownPreprocess() {
        List<Peak> intensPeaks = getIntensPeaks();
        for (Peak p : intensPeaks) {
            int intM2z = (int)(p.getM2z()*10 + 0.5);
            int indexShift = maxShift;
            double sqrtIntens = Math.sqrt(p.getIntensity());
            
            while (indexShift >= 0) {
                double shiftProduct = indexShift*indexShift;
                //double maxShiftProduct = maxShift*maxShift;
                double tempIntens = sqrtIntens* Math.exp(-shiftProduct/2/maxShiftProduct);
                shift(indexShift, intM2z, tempIntens);
                indexShift--;
            }
        }

        // entropy and BUILD_EXP_MASS goes here 
    }

    public int getFinalNumPeaks() {
        return finalNumPeaks;
    }
    protected void xcorrPreprocess() {
        //Iterator <Peak> peakIt = peaks.getPeaks();
        //Iterator <Peak> peakIt = getIntensPeaks(params.getPeakRankThreshold()).iterator();
        //int counter = 0;
        if(boolMasses == null) boolMasses = new boolean[(int) (prcMass+params.getMaxMassShift()+20.5)*PRECISIONFACTOR];
        int myNumPeaks = 0;
        //while (peakIt.hasNext()) {
        List<Peak> peaks = getIntensPeaks(params.getPeakRankThreshold());
        //printPeaksToLog(peaks);
        for(Peak p : peaks) {

            //Peak p = peakIt.next();
            //int intM2z = (int) (10.f*p.getM2z()*MassSpecConstants.DBINWIDTH+0.5f);
            double m2z = p.getM2z();
            if(m2z < prcMass) {
                int intM2z = (int) (PRECISIONFACTOR*m2z + 0.5f);
                double sqrtIntens = Math.sqrt(p.getIntensity());
            //if (sqrtIntens > globMax) {
                int indexShift = maxShift;
                myNumPeaks++;
                while (indexShift >= 0) {
                    double shiftProduct = indexShift*indexShift;
                    //double maxShiftProduct = maxShift*maxShift;
                    // transformedInt as ftemp in PREPROCESS in pep_prob
                    double transformedInt = (sqrtIntens* Math.exp(-shiftProduct/2/maxShiftProduct));
                    int tempInt = (int)(100f*transformedInt/maxTransformedInten + 0.5f);
                    shift(indexShift, intM2z, tempInt); 
                    indexShift--;
                }
            }
           // }
           // counter++;
        }           
        //System.out.print("FNumPeaks: " + myNumPeaks + "\t");
        finalNumPeaks = myNumPeaks;
        //int firstTrue = 0;
        //int lastTrue = 0;
        numTrues =0;
        for(int i = 0; i < boolMasses.length; i++) {
            if (boolMasses[i]) {
                //System.out.println(i);
                if(firstTrue == 0) {
                    firstTrue = i;
                }
                numTrues++;
                lastTrue = i;
            }
        }
        //System.out.print("NumTrues: " + numTrues + "\tFirstTrue: " + firstTrue + "\tLastTrue: " + lastTrue);
        numFragBins = (lastTrue - firstTrue);
        pTrue = (0.0+numTrues)/(lastTrue-firstTrue);
        // set the working aamasses here
        //System.out.println("\tchargeState: " + chargeState + "\tpTrue: " + pTrue + "\t");
    }

    private void writeLogFile() {
        String fileNameBool = "/var/tmp/boolMasses.txt";
        String fileNameIntensity = "/var/tmp/IntensityVals.txt";

        try {

            File file = new File(fileNameBool);
            if (file != null && file.exists()) {
                file.delete();
            }
            file = new File(fileNameIntensity);
            if (file != null && file.exists()) {
                file.delete();
            }
            FileOutputStream outBool = new FileOutputStream(fileNameBool);
            FileOutputStream outIntensity = new FileOutputStream(fileNameIntensity);
            BufferedWriter writerBool = new BufferedWriter(new OutputStreamWriter(outBool));
            BufferedWriter writerIntensity = new BufferedWriter(new OutputStreamWriter(outIntensity));
/*
            String bm = Arrays.toString(boolMasses);
            writerBool.write(bm.replace(",", " ").replace("[", "").replace("]", ""));
            String iv = Arrays.toString(intensityVal);
            writerIntensity.write(iv.replace(",", " ").replace("[", "").replace("]", ""));
*/
            for(boolean b: boolMasses)
            {
                writerBool.write(b+"\t");
                writerBool.newLine();
            }
            for(double d: intensityVal)
            {
                writerIntensity.write(d+"\t");
                writerIntensity.newLine();
            }


            writerBool.close();
            writerIntensity.close();

            outBool.close();
            outIntensity.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    public void printPeaksToLog(float [][] peakarray)
    {
        String location = "/var/tmp/peaksLog.txt";
        try
        {
            BufferedWriter bw = new BufferedWriter(new FileWriter(location));
            for(int i=0; i<peakarray[0].length; i++)
            {
                bw.write(peakarray[0][i]+"\t"+peakarray[1][i]);
                bw.newLine();
            }
            bw.close();
        }
        catch(IOException ioe)
        {

        }
        System.out.println("Peaks log printed at "+location);
    }
    public void printPeaksToLog(List<Peak> peaks)
    {
        String location = "/var/tmp/peaksLog.txt";
        try
        {
            BufferedWriter bw = new BufferedWriter(new FileWriter(location));
            for(Peak peak: peaks)
            {
                bw.write(peak.getM2z()+"\t"+peak.getIntensity());
                bw.newLine();
            }
            bw.close();
        }
        catch(IOException ioe)
        {

        }
        System.out.println("Peaks log printed at "+location);
    }

    public int getFirstTrue() {
        return firstTrue;
    }
    public int getLastTrue() {
        return lastTrue;
    }
    public boolean [] getBoolMasses() {
        if(boolMasses==null) xcorrPreprocess();
        return boolMasses;
    }
    public int getNumTrues() {
        return numTrues;
    }
    public int getNumFragBins() {
        return numFragBins;
    }
    protected void  defaultPreprocess() {
        //int counter = 0;
        for (Peak p : peaks) {
            // intM2z as l in PREPROCESS_SPECTRUM in pep_prob
            int intM2z = (int)(10*(p.getM2z()+0.05f));
            double sqrtIntens = Math.sqrt(p.getIntensity());
            // boolMasses equavelenty to exp_intensity in pep_prob
            int indexShift = maxShift;
            while (indexShift >= 0) {
                shift(indexShift, intM2z, sqrtIntens); 
                indexShift--;
            }
           // counter++;
        }
    }

    private List<Peak> getIntensPeaks(int num) {
        List<Peak> sortedPeaks = peakList.getSortedPeaks(PeakList.SORTBYINTENSITY);
        List<Peak> intensPeaks = new LinkedList<Peak>();
        int numPeaks = sortedPeaks.size();
        int numIntensPeaks = 0;
        while(numIntensPeaks < num && numIntensPeaks < numPeaks) {
            numIntensPeaks++;
            Peak p = sortedPeaks.get(numPeaks-numIntensPeaks);
            intensPeaks.add(p);
        }
        return intensPeaks;
    }
    // TOPDOWN in pep_probe
    protected List<Peak> getIntensPeaks() {
        
        List<Peak> sortedPeaks = peakList.getSortedPeaks(PeakList.SORTBYINTENSITY);
        ArrayList<Peak> intensPeaks = new ArrayList<Peak>(PeakList.DEFAULTNUMPEAKS);
        int numPeaks = sortedPeaks.size();
        double totalIntensity = 0;
        for (Peak p : sortedPeaks) {
            totalIntensity += p.getIntensity();
        }
        double prevAvg = totalIntensity/numPeaks;
        int numIntensPeaks = 0; 
        int numWeakPeaks = 0; // for continues num peaks low than threhhold
        double threshold = 0.025;
         
        while(numIntensPeaks < numPeaks) {
            numIntensPeaks++;
            Peak p = sortedPeaks.get(numPeaks-numIntensPeaks);
            
            intensPeaks.add(p);
            totalIntensity -= p.getIntensity(); 
            double currAvg = totalIntensity/(numPeaks-numIntensPeaks);
            if ((prevAvg - currAvg) <= threshold*prevAvg) {
                if (++numWeakPeaks == 10) {
                   break;
                }
            } else {
                numWeakPeaks = 0;
            }
            prevAvg = currAvg; 
        }
        return intensPeaks;
    }
    public double getPTrue() {
        return pTrue;
    }


    public float[] getMassCorrFloat()
    {
        float[] result = new float[massCorr.length];
        for(int i=0; i<massCorr.length; i++)
        {
            result[i] = (float)massCorr[i];
        }
        return result;
    }
    public float[] getMassSumFloat()
    {
        float[] result = new float[massSum.length];
        for(int i=0; i<result.length;i++)
        {
            result[i] = (float) massSum[i];
        }
        return result;
    }
    public void dumpBoolMass()
    {
        boolMasses = null;
    }
}



