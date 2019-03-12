/**
 * @file ProcessedPeakList.java
 * This is the source file for edu.scripps.pms.mspid.ProcessedPeakList
 * @author Tao Xu
 * @date $Date: 2009/09/01 05:24:59 $
 */
package edu.scripps.pms.mspid;

import edu.scripps.dia.LibrarySpectra;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.spectrum.*;
import gnu.trove.*;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.regression.SimpleRegression;

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
    private static int ID_SEED=0;
    private int id ;

    protected PeakList peakList;
    protected ArrayList<Peak> peaks;
    protected SearchParams params;
    //protected double [] intensityVal = new double[MAXNUMPEAKS];
    protected double [] intensityVal; // = new double[MAXNUMPEAKS];
    protected double [] intensityVal2;

    // as exp_intensity in pep_prob
    //protected boolean [] boolMasses = new boolean[MAXNUMPEAKS]; // boolean representation of masses
    protected boolean [] boolMasses; // = new boolean[MAXNUMPEAKS]; // boolean representation of masses
    protected TIntHashSet massSet;
    protected TIntIntHashMap massPeakIdMap;

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
    private List<Peak> topPeaks ;
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


    public ProcessedPeakList(PeakList peaklist, Zline z, SearchParams sp, MassCalculator mc) {
        init(peaklist,z,sp,mc);
    }

    private void init(PeakList peaklist, Zline z, SearchParams sp, MassCalculator mc)
    {
        id = ID_SEED++;
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

    public TIntIntHashMap getMassPeakIdMap() {
        if(massPeakIdMap ==null) xcorrPreprocess();
        return massPeakIdMap;
    }
    public ScoredPeptideHit prelimScoreCorrelation(ProcessedPeakList peakList)
    {
        int numTheroticPeaks =0;
        int numPeaksMatched =0;
        int i=0;
        SimpleRegression sr = new SimpleRegression();
        // boolean [] boolMass = peakList.getBoolMasses();
        // boolean [] mybools = this.getBoolMasses();

        //TIntHashSet altSet = peakList.getMassSet();
        //TIntHashSet mySet = this.getMassSet();

        //List<Peak> peaks = getIntensPeaks(params.getPeakRankThreshold());
        TIntIntHashMap myMap = getMassPeakIdMap();
        List<Peak> altpeaks = peakList.getIntensPeaks(params.getPeakRankThreshold());
        numTheroticPeaks = altpeaks.size();

        for(Peak p: altpeaks)
        {
            int mass  = (int)(p.getM2z() *PRECISIONFACTOR+ 0.5f);
            if(myMap.containsKey(mass))
            {
                int myPeakID = myMap.get(mass);
                //  System.out.println("+++"+altPeakID);
                double intensity =getIntensityFromTopPeak(myPeakID);
                sr.addData(p.getIntensity(),intensity);
                numPeaksMatched++;
            }
        }



    /*    TIntIntHashMap altMap = peakList.getMassPeakIdMap();
        TIntIntHashMap myMap = getMassPeakIdMap();
        TIntHashSet idSet = new TIntHashSet();

        numTheroticPeaks = altMap.size();
        for(TIntIntIterator titr = myMap.iterator(); titr.hasNext(); )
        {
            titr.advance();
            int mass = titr.key();
            int id = titr.value();
            if(altMap.containsKey(mass))
            {
                numPeaksMatched++;
                if(!idSet.contains(id))
                {
                    double myIntensity = getIntensityFromTopPeak(id);
                    int altID = altMap.get(mass);
                    double altIntenisty = peakList.getIntensityFromTopPeak(altID);
                    sr.addData(myIntensity,altIntenisty);
                    idSet.add(id);
                }
            }

        }*/


     /*   numTheroticPeaks = peakList.getTopPeaksSize();
        for(Peak p: peaks)
        {
            int mass  = (int)(p.getM2z() *PRECISIONFACTOR+ 0.5f);
            if(altMap.containsKey(mass))
            {
                int altPeakID = altMap.get(mass);
              //  System.out.println("+++"+altPeakID);
                double intensity =peakList.getIntensityFromTopPeak(altPeakID);
                sr.addData(p.getIntensity(),intensity);
                numPeaksMatched++;
            }
        }*/

        /*for(int j=firstTrue; j<lastTrue; j++)
        {
            if(boolMass[j] )
            {
                numTheroticPeaks++;
              if(mybools[j]) numPeaksMatched++;
            }
        }*/
        //peakList.dumpBoolMass();
        //this.dumpBoolMass();
        double probability = DistributionCalculator.getBinomialSum(getPTrue(), numTheroticPeaks, numPeaksMatched);
        ScoredPeptideHit sph = new ScoredPeptideHit(probability);
        sph.setPrcMass(peakList.prcMass);
        sph.setrSqaured(sr.getRSquare());
        sph.setRetTime(peakList.getPeakList().getRetentionTime());
        return sph;
    }


    public ScoredPeptideHit prelimScoreCorrelation(LibrarySpectra spectra)
    {
        int numTheroticPeaks =0;
        int numPeaksMatched =0;
        int i=0;
        SimpleRegression sr = new SimpleRegression();

        TIntIntHashMap myMap = getMassPeakIdMap();
        List<Peak> altpeaks = spectra.getIntensPeaks(params.getPeakRankThreshold());
        numTheroticPeaks = altpeaks.size();

        for(Peak p: altpeaks)
        {
            int mass  = (int)(p.getM2z() *PRECISIONFACTOR+ 0.5f);
            if(myMap.containsKey(mass))
            {
                int myPeakID = myMap.get(mass);
                //  System.out.println("+++"+altPeakID);
                double intensity =getIntensityFromTopPeak(myPeakID);
                sr.addData(p.getIntensity(),intensity);
                numPeaksMatched++;
            }
        }

        double probability = DistributionCalculator.getBinomialSum(getPTrue(), numTheroticPeaks, numPeaksMatched);
        ScoredPeptideHit sph = new ScoredPeptideHit(probability);
        sph.setNumMatchedPeaks(numPeaksMatched);
        sph.setMs2CompareValues(spectra.scan,spectra.scan,spectra.filename);
        sph.setPrcMass(spectra.mz);
        sph.setrSqaured(sr.getRSquare());
        sph.setNumPeaks(altpeaks.size());
        sph.setLibrarySpectra(spectra);
        sph.setRetTime(spectra.retTime);
        return sph;
    }

    public ScoredPeptideHit dotProduct(LibrarySpectra spectra)
    {
        double score = calculateDotProductDerived(spectra, 1, 0,false);

        ScoredPeptideHit sph = new ScoredPeptideHit(score);
        sph.setMs2CompareValues(spectra.scan,spectra.scan,spectra.filename);
        sph.setPrcMass(spectra.mz);
        sph.setNumPeaks(spectra.getMzList().size());
        sph.setLibrarySpectra(spectra);
        sph.setRetTime(spectra.retTime);
        sph.setIsDecoy(spectra.isDecoy());
        return sph;
    }


    public ScoredPeptideHit normalizedDotProduct(LibrarySpectra spectra, float x , float y)
    {
        double score = calculateDotProductDerived(spectra, x, y,true);

        ScoredPeptideHit sph = new ScoredPeptideHit(score);
        sph.setMs2CompareValues(spectra.scan,spectra.scan,spectra.filename);
        sph.setPrcMass(spectra.mz);
        sph.setNumPeaks(spectra.getMzList().size());
        sph.setLibrarySpectra(spectra);
        sph.setRetTime(spectra.retTime);
        sph.setIsDecoy(spectra.isDecoy());
        return sph;
    }

    public ScoredPeptideHit pearsonsCorrelation(LibrarySpectra spectra)
    {
        double score = calcPearsonsCorrelation(spectra);

        ScoredPeptideHit sph = new ScoredPeptideHit(score);
        sph.setMs2CompareValues(spectra.scan,spectra.scan,spectra.filename);
        sph.setPrcMass(spectra.mz);
        sph.setNumPeaks(spectra.getMzList().size());
        sph.setLibrarySpectra(spectra);
        sph.setRetTime(spectra.retTime);
        sph.setIsDecoy(spectra.isDecoy());
       // System.out.println(">> SCORE "+score);
        return sph;
    }

    public ScoredPeptideHit spearmansCorrelation(LibrarySpectra spectra)
    {
        double score = calcSpearmansCorrelation(spectra);

        ScoredPeptideHit sph = new ScoredPeptideHit(score);
        sph.setMs2CompareValues(spectra.scan,spectra.scan,spectra.filename);
        sph.setPrcMass(spectra.mz);
        sph.setNumPeaks(spectra.getMzList().size());
        sph.setLibrarySpectra(spectra);
        sph.setRetTime(spectra.retTime);
        sph.setIsDecoy(spectra.isDecoy());
        return sph;
    }

    public double calcPearsonsCorrelation(LibrarySpectra spectra)
    {
        PeakList libpeaks = spectra.getPeakList();
        float min = libpeaks.getMinM2z()<peakList.getMinM2z()? (float)libpeaks.getMinM2z(): (float)peakList.getMinM2z();
        float max = libpeaks.getMaxM2z()>peakList.getMaxM2z()? (float)libpeaks.getMaxM2z(): (float)peakList.getMaxM2z();
        PearsonsCorrelation correlation = new PearsonsCorrelation();
        float[] yArray = spectra.getPeakList().generateBinSpectra(0,(float)fragTolerance,
                0,min,max).toNativeArray();
        float[] xArray = this.peakList.generateBinSpectra(0,(float)fragTolerance,
                0,min,max).toNativeArray();
        double[] xArrayd = new double[xArray.length];
        double[] yArrayd = new double[yArray.length];
        for(int i=0; i<xArray.length; i++)
        {
            xArrayd[i] = xArray[i];
        }
        for(int i=0; i<yArray.length; i++)
        {
            yArrayd[i] = yArray[i];
        }
        return correlation.correlation(xArrayd,yArrayd);
    }

    public double calcSpearmansCorrelation(LibrarySpectra spectra)
    {
        PeakList libpeaks = spectra.getPeakList();
        float min = libpeaks.getMinM2z()<peakList.getMinM2z()? (float)libpeaks.getMinM2z(): (float)peakList.getMinM2z();
        float max = libpeaks.getMaxM2z()>peakList.getMaxM2z()? (float)libpeaks.getMaxM2z(): (float)peakList.getMaxM2z();
        SpearmansCorrelation correlation = new SpearmansCorrelation();
        float[] yArray = spectra.getPeakList().generateBinSpectra(0,(float)fragTolerance,
                0,min,max).toNativeArray();
        float[] xArray = this.peakList.generateBinSpectra(0,(float)fragTolerance,
                0,min,max).toNativeArray();
        double[] xArrayd = new double[xArray.length];
        double[] yArrayd = new double[yArray.length];
        for(int i=0; i<xArray.length; i++)
        {
            xArrayd[i] = xArray[i];
        }
        for(int i=0; i<yArray.length; i++)
        {
            yArrayd[i] = yArray[i];
        }
        return correlation.correlation(xArrayd,yArrayd);
    }



    public double calculateDotProductDerived(LibrarySpectra spectra,
                                             float x, float y, boolean isNormalized) {
        // TODO: Make sure if any of binSpectra is empty, score i NaN! and return!
        // now, calculate dot product score
        PeakList libpeaks = spectra.createPeakListCopy();
        float min = libpeaks.getMinM2z()<peakList.getMinM2z()? (float)libpeaks.getMinM2z(): (float)peakList.getMinM2z();
        float max = libpeaks.getMaxM2z()>peakList.getMaxM2z()? (float)libpeaks.getMaxM2z(): (float)peakList.getMaxM2z();
       /* if(spectra.isDecoy())
        {
            System.out.print(">> ");
            for(int i=0; i<spectra.getMzList().size(); i++)
            {
                System.out.print(spectra.getMzList().get(i)+", "+spectra.getIntensityList().get(i)+" ");
            }
            System.out.println();
        }*/
        float[] yArray =libpeaks.generateBinSpectra(0,(float)fragTolerance,
                0,min,max).toNativeArray();
        float[] xArray = this.peakList.generateBinSpectra(0,(float)fragTolerance,
                0,min,max).toNativeArray();

        float dot_product_alpha_beta = 0,
                dot_product_alpha_alpha = 0,
                dot_product_beta_beta = 0;
        float binSize = (float)(2 * fragTolerance),
                mz = (float)peakList.getMinM2z() + binSize;
        for (int i = 0; i < xArray.length; i++) {
            float mz_peakA = mz,
                    mz_peakB = mz,
                    intensity_peakA = xArray[i],
                    intensity_peakB = yArray[i];
            boolean control = false;
            if (intensity_peakA == 0 && intensity_peakB == 0) {
                control = true;
            }
            if (!control) {
                if (i == xArray.length - 1) {// make sure last one is not empty..
                    mz_peakA = mz_peakA + 0.0000001f;
                    mz_peakB = mz_peakB + 0.0000001f;
                }
                float alpha =(float)( Math.pow(intensity_peakA, x) * Math.pow(mz_peakA, y)),
                        beta =(float)( Math.pow(intensity_peakB, x) * Math.pow(mz_peakB, y));
                dot_product_alpha_beta = dot_product_alpha_beta + (float) (alpha * beta);
                dot_product_alpha_alpha = dot_product_alpha_alpha + (float) (alpha * alpha);
                dot_product_beta_beta = dot_product_beta_beta + (float) (beta * beta);
                mz = mz + binSize;
            }
        }
        float score;
        if (isNormalized) {
            float normalized_dot_product = (float) dot_product_alpha_beta / (float) (Math.sqrt(dot_product_alpha_alpha * dot_product_beta_beta));
            score = normalized_dot_product;
        } else {
            score = dot_product_alpha_beta;
        }
        return score;
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


    public void preprocess(int mode) {
        TimeUtils timer = new TimeUtils();
        //timer.startTiming();
        switch (mode) {

            case XCORRPREPROCESS: xcorrPreprocess();
               // writeLogFile();
            break;

        }
        //System.out.print("Scan#: " + peakList.getLoscan() + "\t+" + zline.getChargeState() + "\tNumPeaks: " + peakList.numPeaks()+ "\tm2z: " + zline.getM2z() + "\tptrue: " + pTrue);
    }

    protected void shift(int indexShift, int intM2z, int id) {

        int leftIndex = intM2z - indexShift;
        int rightIndex = intM2z + indexShift;
        //boolMasses[leftIndex] = true;
        //boolMasses[rightIndex] = true;
        ///massSet.add(leftIndex);
        massPeakIdMap.put(leftIndex,id);
        massPeakIdMap.put(rightIndex,id);

        //massSet.add(rightIndex);
        if(leftIndex<firstTrue)
        {
            firstTrue = leftIndex;
        }
        if(rightIndex>lastTrue)
        {
            lastTrue = rightIndex;
        }

/*        if (intensityVal[leftIndex] <= transformedIntensity) {
            intensityVal[leftIndex] = transformedIntensity;
        }

        if (intensityVal[rightIndex] <= transformedIntensity) {
            intensityVal[rightIndex] = transformedIntensity;
        }*/
    }



    protected void shift(int indexShift, int intM2z, double transformedIntensity) {
        
        int leftIndex = intM2z - indexShift;
        int rightIndex = intM2z + indexShift;
        //boolMasses[leftIndex] = true;
        //boolMasses[rightIndex] = true;
        ///massSet.add(leftIndex);
     //   massPeakIdMap.put(leftIndex,(float)transformedIntensity);
     //   massPeakIdMap.put(rightIndex,(float)transformedIntensity);

        //massSet.add(rightIndex);
        if(leftIndex<firstTrue)
        {
            firstTrue = leftIndex;
        }
        if(rightIndex>lastTrue)
        {
            lastTrue = rightIndex;
        }

/*        if (intensityVal[leftIndex] <= transformedIntensity) {
            intensityVal[leftIndex] = transformedIntensity;
        } 
           
        if (intensityVal[rightIndex] <= transformedIntensity) {
            intensityVal[rightIndex] = transformedIntensity;
        }*/
    }


    public int getFinalNumPeaks() {
        return finalNumPeaks;
    }
    protected void xcorrPreprocess() {
        //Iterator <Peak> peakIt = peaks.getPeaks();
        //Iterator <Peak> peakIt = getIntensPeaks(params.getPeakRankThreshold()).iterator();
      //  if(boolMasses == null) boolMasses = new boolean[(int) (prcMass+params.getMaxMassShift()+20.5)*PRECISIONFACTOR];
        //if(massSet == null) massSet = new TIntHashSet();
         //massSet = new TIntHashSet();
        massPeakIdMap = new TIntIntHashMap();
        int myNumPeaks = 0;
        //while (peakIt.hasNext()) {
        if(topPeaks==null)
         topPeaks = getIntensPeaks(params.getPeakRankThreshold());
        //printPeaksToLog(peaks);
        firstTrue = Integer.MAX_VALUE;
        lastTrue = Integer.MIN_VALUE;
        int peakID = 0;

        for(Peak p : topPeaks) {

            //Peak p = peakIt.next();
            //int intM2z = (int) (10.f*p.getM2z()*MassSpecConstants.DBINWIDTH+0.5f);
            double m2z = p.getM2z();
            if(m2z < prcMass) {
                int intM2z = (int) (PRECISIONFACTOR*m2z + 0.5f);
            //if (sqrtIntens > globMax) {
                int indexShift = maxShift;
                myNumPeaks++;
                while (indexShift >= 0) {
                    //double maxShiftProduct = maxShift*maxShift;
                    // transformedInt as ftemp in PREPROCESS in pep_prob
                    shift(indexShift, intM2z,peakID );

                    //shift(indexShift, intM2z, p.getIntensity());
                    indexShift--;
                }

            }
            peakID++;
           // }
        }           
        //System.out.print("FNumPeaks: " + myNumPeaks + "\t");
        finalNumPeaks = myNumPeaks;
        //int firstTrue = 0;
        //int lastTrue = 0;
       // numTrues =0;
        //firstTrue =0;
   /*     for(int i = 0; i < boolMasses.length; i++) {
            if (boolMasses[i]) {
                //System.out.println(i);
                if(firstTrue == 0) {
                    firstTrue = i;
                }
                numTrues++;
                lastTrue = i;
            }
        }*/
        numTrues = massPeakIdMap.size();
        //System.out.print("NumTrues: " + numTrues + "\tFirstTrue: " + firstTrue + "\tLastTrue: " + lastTrue);
        numFragBins = (lastTrue - firstTrue);
        pTrue = (0.0+numTrues)/(lastTrue-firstTrue);
        // set the working aamasses here
        //System.out.println("\tchargeState: " + chargeState + "\tpTrue: " + pTrue + "\t");
    }
    public double getIntensityFromTopPeak(int peakId)
    {
        return topPeaks.get(peakId).getIntensity();
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

    public TIntHashSet getMassSet()
    {
        if(massSet == null) xcorrPreprocess();
        return massSet;
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
        massSet = null;
        massPeakIdMap = null;
        topPeaks = null;
    }
    public int getID()
    {
        return id;
    }

    public int getTopPeaksSize()
    {
        return topPeaks.size();
    }

}



