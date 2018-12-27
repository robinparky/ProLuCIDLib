/**
 * @file SearchResult.java
 * This is the source file for edu.scripps.pms.mspid.SearchResult
 * @author Tao Xu
 * @date $Date: 2015/04/30 21:01:12 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.mspid.ProcessedPeakList;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.stat.StatCalc;

import java.text.DecimalFormat;
import java.util.*;

public class SearchResult {
    /* for binomial and poisson
    private final int NUM = 36;
    private int [] freq = new int[NUM];
    private int [] freqall = new int[NUM];
    */
    public static int NUMSCORED;
    private static int NUMFINALRESULT = 5;
    private static final char DELIMITER = '\t';
    private double retTime;
    //private LinkedList<PeptideHit> hits = new LinkedList<PeptideHit>(); 
    private ArrayList<PeptideHit> hits;
    private List<ScoredPeptideHit> scoredHits;
    private List<ScoredPeptideHit> finalResult = null;
    private boolean isSorted = true;
//    private int totalPeptideLength = 0;
    private int numPeaks = 0;  // total number of peaks from all matched peptides
    private int totalPeptideLength = 0;
    private int numPeaksMatched = 0;
    private final int minNumPeaksMatched;
    private final int minPeptideLength;
    private int chargeState;
    private String ms2Name;
    //private final int effectChargeState; // chargeState of ppl
    private ProcessedPeakList ppl;
    private SearchParams params;
    private int numPeptidesMatched = 0;
    private MassCalculator mc;
    private long searchTime;
    private String hostName;
    private PeakList peaks;
    private Protease protease;
    private int enzymeSpecificity;
    public static DecimalFormat threeDigits = new DecimalFormat("0.000");
    private static DecimalFormat fourDigits = new DecimalFormat("0.0000");
    private static DecimalFormat deltaCNDigits = new DecimalFormat("0.0000E0");

    private static DecimalFormat fiveDigits = new DecimalFormat("0.00000");
    private static DecimalFormat twoDigits = new DecimalFormat("0.00");
    private ScoredPeptideHit bestSecondaryScore;
    private boolean sortByXcorr = false;
    private int primaryScoreType;     
    private int secondaryScoreType;     
    private double primaryScoreDeviation = 1;
    private double primaryScoreMean = 0;
//    private static DistributionCalculator dc = new DistributionCalculator();

    public SearchResult(ProcessedPeakList ppl) {
        this.ppl = ppl;
        params = ppl.getSearchParams();
        protease = params.getProtease();
        
        enzymeSpecificity = params.getEnzymeSpecificity();
        chargeState  = ppl.getZline().getChargeState();       
        //effectChargeState  = ppl.getZline().getChargeState() > 2? 3 : 2;       
        peaks = ppl.getPeakList();
        mc = ppl.getMassCalculator();
        minNumPeaksMatched = params.getMinMatch();
        minPeptideLength = ppl.getSearchParams().getMinimumPeptideLength() - 2;
        NUMSCORED = ppl.getSearchParams().getCandidatePeptideThreshold();
        hits = new ArrayList<PeptideHit>(NUMSCORED+1);
        scoredHits = new ArrayList<ScoredPeptideHit>(NUMSCORED);
        primaryScoreType = ppl.getSearchParams().getPrimaryScoreType(); 
        secondaryScoreType = ppl.getSearchParams().getSecondaryScoreType(); 
        sortByXcorr = (primaryScoreType != 0); // ie, primaryScoreType not probability score
    }
    public void setHostName(String host) {
        hostName = host;
    }
    public static void setNumOutput(int num) {
        NUMFINALRESULT = num;
    } 
    public void setSearchTime(long timeUsedInMilliSeconds) {
        searchTime = timeUsedInMilliSeconds;
    }
    protected void ignoreModifiedPeptideHit() {
        ScoredPeptideHit  topModified = null;
        ScoredPeptideHit topUnModified = null;
        for(Iterator<ScoredPeptideHit> it = scoredHits.iterator(); it.hasNext();) {
            ScoredPeptideHit hit = it.next();
            if(hit.isModified()) {
                if(topModified == null) {
                    topModified = hit;
                }
            } else {
                if(topUnModified == null) {
                    topUnModified = hit;
                }            
            }
            if(topModified != null && topUnModified != null) {
                break;
            }
        } 
        if(topModified != null && topUnModified != null) {
            //if(topModified.isIgnorable(topUnModified)) {}
            // need to implement the criteria to determine ignore or not
        } 
    } 
    public List<ScoredPeptideHit> getFinalResult() {
       

        ArrayList<ScoredPeptideHit> finalList = new ArrayList<ScoredPeptideHit>();
        if(scoredHits.size() == 0) {
            return finalList;
        }
        String topseq = scoredHits.get(0).getOriginalSequence();

        int numAdded = 0;
        
        for(int i = 0; i < scoredHits.size(); i++) {
            ScoredPeptideHit sph = scoredHits.get(i);
            finalList.add(scoredHits.get(i));
            numAdded++;
//System.out.println(topseq + "\t" + sph.getOriginalSequence() + "\t" + NUMFINALRESULT + "\t" + sph.getPrimaryRank() + "\t" + !topseq.equals(sph.getOriginalSequence()));
//System.out.println(NUMFINALRESULT + "\t" + i + "\t" + (i < NUMFINALRESULT) + "\t" +  !topseq.equals(sph.getOriginalSequence()));
            if(numAdded >= NUMFINALRESULT && (!topseq.equals(sph.getOriginalSequence()))) {
                break; 
            }
        }
        //finalResult.add(bestSecondaryScore);
        return finalList;

        // the following was trying to do something special with modifcation hits
        /*
        ArrayList<ScoredPeptideHit> finalList = new ArrayList<ScoredPeptideHit>();
        for(int i = 0; finalList.size()< NUMFINALRESULT && i < scoredHits.size(); i++) {
            ScoredPeptideHit hit = scoredHits.get(i);
            if(!ignoreModifiedPeptideHit() || !hit.isModified()) {
                finalList.add(scoredHits.get(i));
            } 
        }
        //finalResult.add(bestSecondaryScore);
        if(finalList.size() > 0) {
            finalList.get(0).setPrimaryRank(1);
        }
        return finalList;
        */
    }
    public double getDeltaCn(ScoredPeptideHit sph) {
        /*
        if(tscore) {
            // use t score to replace deltaCN
            return (sph.getPrimaryScore() - primaryScoreMean)/primaryScoreDeviation;
        } else {
            double bestScore = finalResult.get(0).getPrimaryScore();       
            if(bestScore == 0) {
                return 0;
            }
            return (bestScore - sph.getPrimaryScore())/bestScore;
        }
        */
        double bestScore = finalResult.get(0).getPrimaryScore();       
        if(bestScore == 0) {
            return 0;
        }
        return (bestScore - sph.getPrimaryScore())/bestScore;

    } 
    public ScoredPeptideHit getTopHit() {
        if(finalResult != null && finalResult.size() > 0) {
            return finalResult.get(0);
        } else {
            return null;
        }
    }
    public int getNumPeptidesMatched() {
        return numPeptidesMatched;
    }

    public ProcessedPeakList getProcessedPeakList() {
        return ppl;
    }
    public String outputResults() {
        //calcScores();
        
        StringBuffer result = new StringBuffer(8000);
        appendSline(result);
        appendMlines(result); // L lines are also appended
        result.append('\n');
/*
System.out.println("numPeaks: " + numPeaksMatched + "\tnumPeptides: " + numPeptidesMatched);
double lamda = numPeaksMatched/(numPeptidesMatched+0.0); 
System.out.println("lamda: \t" + lamda);
System.out.println("numPeaksMatched\tlength " + NUM/2 + "\tallLength\t");

for(int i = 0; i < freq.length; i++) {
    System.out.println(i + "\t" + freq[i] + "\t" + freqall[i]);
}
*/
        return result.toString();
        
    }
    

    private void appendMlines(StringBuffer result) {
        if(finalResult == null || finalResult.size() < 1) { return; }
        for(ScoredPeptideHit p : finalResult) {
            appendMline(result, p);
        }
    }
    public static final String MLINEHEADER = "MLINE\tPSCORE\tSCANHI\tSCANLOW\tPRCMASS\tR-SQUARED\tRET-TIME\tIS-DECOY\tTARGET-PATH";

    private void appendMLineSpectraComparison(StringBuffer result, ScoredPeptideHit sph)
    {
        result.append("M");
        result.append(DELIMITER);
        result.append(sph.getPScore());
        result.append(DELIMITER);
        result.append(sph.getScanHi());
        result.append(DELIMITER);
        result.append(sph.getScanLow());
        result.append(DELIMITER);
        result.append(sph.getPrcMass());
        result.append(DELIMITER);
        result.append(sph.getrSqaured());
        result.append(DELIMITER);
        result.append(sph.getRetTime());
        result.append(DELIMITER);
        String name = sph.isDecoy()? ("Decoy\t"+sph.getMs2Filename()) :("\t"+sph.getMs2Filename());
        result.append(name);
        result.append("\n"); // Validation state



    }

    private void appendMlineIndexedSpectra(StringBuffer result, ScoredPeptideHit sph) {
        // Format of M line:
        // "M\tRankByXcorr\tRankBySp\tCalculatedMass\tDeltaCN\tXcorr\t
        // Sp\tNumMatchedIons\tNumExpectedIons\tPeptideSequence\tValidateState\n"
        if (sph == null) return;
        result.append("M");
        result.append(DELIMITER);
        result.append(sph.getPrimaryRank());
        result.append(DELIMITER);
        result.append(sph.getPrimaryRank());
        result.append(DELIMITER);
        result.append(fiveDigits.format(sph.getPrcMass()));
        result.append(DELIMITER);
        result.append(fiveDigits.format(getDeltaCn(sph))); //deltaCn
        result.append(DELIMITER);
        result.append(fourDigits.format(sph.getPrimaryScore())); // XCorr
        result.append(DELIMITER);
        //result.append(threeDigits.format(sph.getSecondaryScore())); // XCorr
        result.append(sph.getSecondaryScore()); // for probablity scores, cannot use format
        result.append(DELIMITER);
      //  result.append(sph.getNumPeaksMatched()); //numPeaksMatched
        result.append(DELIMITER);
        result.append(sph.getNumPeaksMatched()); //numPeaks
        result.append(DELIMITER);
        result.append(sph.getNumPeaks()); //numPeaks
        result.append(DELIMITER);
        result.append(sph.getExtendedSequence()); // Peptide sequence
        result.append(DELIMITER);
     //   result.append(sph.getLibrarySpectra().filename);
      //  result.append(DELIMITER);
/*
        result.append(ppl.getZline().getScanInfo());
        result.append(DELIMITER);*/
        result.append("U\n"); // Validation state

        int locusType = params.getLocusType();
        //append Llines

        List<String> accessionList = sph.getLibrarySpectra().getAccessionList();
        List<String> descriptionList = sph.getLibrarySpectra().getDescriptionList();

        for(int i=0; i<accessionList.size(); i++)
        {
            result.append("L\t").append(accessionList.get(i)).append("\t")
                    .append(0).append("\t")
                    .append(sph.getSequence()).append("\n");
        }

       /*
        for(PeptideHit p : sph.getPeptideHits()) {



            if(locusType == 1) {
                result.append("L\t" + p.getSequestLikeAccession());
            } else {
                result.append("L\t" + p.getAccession());
            }

            List<String> defList = p.getDefList();
            for(Iterator<String> itr=defList.iterator(); itr.hasNext(); ) {
                String eachDef = itr.next();

                result.append("L\t").append(eachDef).append("\t")
                        .append(p.getStart()).append("\t")
                        .append(p.getExtraExtendedSequence()).append("\n");
            }
        }
        */
    }


    private void appendMline(StringBuffer result, ScoredPeptideHit sph) {
        // Format of M line:
        // "M\tRankByXcorr\tRankBySp\tCalculatedMass\tDeltaCN\tXcorr\t
        // Sp\tNumMatchedIons\tNumExpectedIons\tPeptideSequence\tValidateState\n"

        appendMlineIndexedSpectra(result,sph);
            return;



    }

    public static final String SLINEHEADER = "SLINE\tSCANHI\tSCANLOW\tCHARGE\tHOSTNAME\tPRCMASS\tTOTAL-INTENSITY\tP-TRUE" +
            "\tNUM-SPECTRA-MATCHED\tRET-TIME\tMSNAME";

    private void appendSline(StringBuffer result) {

        // Format of S line:
        // "S\tLowScan\tHighScan\tChargeState\tProcessTime\tSever\t
        // ObservedMass\tTotalIntensity\tLowestSp\tNumSequencesMatched
        result.append("S");
        result.append(DELIMITER);
        result.append(ppl.getPeakList().getHiscan());
        result.append(DELIMITER);
        result.append(ppl.getPeakList().getLoscan());
        result.append(DELIMITER);
        result.append(chargeState);
        result.append(DELIMITER);
        result.append(searchTime);
        result.append(DELIMITER);
        result.append(hostName + "_" + Thread.currentThread().getName());
        result.append(DELIMITER);
        result.append(fiveDigits.format(ppl.getZline().getM2z()));
        result.append(DELIMITER);
        result.append(twoDigits.format(peaks.getTotalIntensity()));
        result.append(DELIMITER);
        result.append(fourDigits.format(ppl.getPTrue())); // for ptrue 
        result.append(DELIMITER);
        result.append(numPeptidesMatched);
        result.append(DELIMITER);
        result.append(ppl.getPeakList().getRetentionTime());
        result.append(DELIMITER);
        //result.append(ms2Name);



        result.append("\n");
    }
    public void calcScores() {
        setPScoreRank(finalResult);
       // calcScores(null);
        }


    private void calcPrimaryScoreDeviation() {
        List<ScoredPeptideHit> resultList = scoredHits;
        //List<ScoredPeptideHit> resultList = finalResult; 
        Iterator<ScoredPeptideHit> it = resultList.iterator();
        int numElement = 0;
        String topseq = null;
        while(it.hasNext()) {
            ScoredPeptideHit p = it.next();
            if(topseq == null) {
                topseq = p.getOriginalSequence();
            }
            //if(p.getPrimaryRank() != 1 && (!topseq.equals(p.getOriginalSequence()))) {
            if(p.getXCorrRank() != 1 && (!topseq.equals(p.getOriginalSequence()))) {
                primaryScoreMean += p.getXCorr();
//System.out.println("PrimaryScore: " + p.getPrimaryScore());
                numElement++;
            }
        }
        if(numElement > 2) { 
            primaryScoreMean /= numElement;
            it = resultList.iterator();
            while(it.hasNext()) {
                ScoredPeptideHit p = it.next();
                if(p.getXCorrRank() > 1 && (!topseq.equals(p.getOriginalSequence()))) {
                    double diff = p.getXCorr() - primaryScoreMean;
                    primaryScoreDeviation  += diff*diff;
                }
            }
            primaryScoreDeviation = Math.sqrt(primaryScoreDeviation/(numElement-1));
            if(primaryScoreDeviation == 0) {
                primaryScoreDeviation = 1;
            }
        } 

//       System.out.println("Number of Elements for zscore: " + numElement);        
       // System.out.println("PrimaryScoreMean: " + primaryScoreMean + "\tPrimaryScoreDeviation: " + primaryScoreDeviation);
            
    }
    private void setPScoreRank(List<ScoredPeptideHit> scoredList) {

        int rank = 0;
        double lastScore = 100;
        for(ScoredPeptideHit p : scoredList) {
            double pScore = p.getPScore();
            if(pScore < lastScore) {
                rank++;
                lastScore = pScore;
            }
            p.setPScoreRank(rank);
        }
    }
    private void setXCorrRank(List<ScoredPeptideHit> scoredList) {

        int rank = 0;
        double lastScore = 100;
        for(ScoredPeptideHit p : scoredList) {
            double xcorr = p.getXCorr();
            if(xcorr < lastScore) {
                rank++;
                lastScore = xcorr; 
            } 
            p.setXCorrRank(rank);
        }

    }
    private void chebyshev() {
        int numPeptideHits = hits.size();
        // j in CHEBYSHEV
        int avgNumPeaksMatched = Math.round(numPeaksMatched/(float)numPeptideHits);
        // k in CHEBYSHEV
        int avgPeptideLength = Math.round(totalPeptideLength/(float)numPeptideHits);
                            
        int avgPeaks = 2*(avgPeptideLength-1)*(chargeState-1);
        float avgHypergeometry = ScoreCalculator.hypergeometry(
                 numPeaks, numPeaksMatched, avgPeaks, avgNumPeaksMatched);
        float avgScore = -(float) Math.log(avgHypergeometry);


        // calculate the best five scores here 
    }

    protected PeptideHit createPeptideHit(Peptide p) {
        if(ppl.isDeCharged()) {
            return new DeChargedPeptideHit(p);
        } else {
           return new PeptideHit(p);
        }
    }
    protected ModifiedPeptideHit createModifiedPeptideHit(Peptide p, Modifications m) {
        if(ppl.isDeCharged()) {
            return new DeChargedModifiedPeptideHit(p, m);
        } else {
           return new ModifiedPeptideHit(p, m);
        }
    }
    public void addPeptideHit(Fasta f, int start, int end) {
        if(isValidHit(f, start, end)) {
            //addFinalPeptideHit(createPeptideHit(new Peptide(f, start, end)));
            addPeptideHit(createPeptideHit(new Peptide(f, start, end)));
        }
    }
    public void addPeptideHit(Fasta f, int start, int end, Modifications m) {
        if(isValidHit(f, start, end)) {
            //addFinalPeptideHit(createModifiedPeptideHit(new Peptide(f, start, end), m));
            addPeptideHit(createModifiedPeptideHit(new Peptide(f, start, end), m));
        }
    }


    public void addPeptideHit(PeptideHit p) {
        //if(isValidHit(p)) {
            addFinalPeptideHit(p);
        //}
    }
    public void addPeptideHit(ModifiedPeptideHit p) {

        //if(isValidHit(p)) {
            // check enzyme specificity 
            for(Iterator<ModifiedPeptideHit> it = p.getAllModifiedPeptideHits(); it.hasNext();) {
           //System.out.println("in add modified peptide hit"); 
                addFinalPeptideHit(it.next());
            }
        //}
    }

    // for semi-blind modification search
    public void addPeptideHit(Peptide p, DiffMod m) {
        
        if(isValidHit(p)) {
            int length = p.getLength();
            for(int i = 0; i < length; i++) {
                ModifiedPeptideHit mph = new ModifiedPeptideHit(p);
               
                mph.setDiffMod(i, m);
                addFinalPeptideHit(mph);
            }
        }
    }

    private boolean isValidHit(Fasta f, int start, int end) {
        if((end-start) > minPeptideLength) {
            // check enzyme specificity 
            if(protease == null || protease.checkEnzymeSpecificity(f, start, end) >= enzymeSpecificity) { 
                return true;
            }
        }
        return false;

    }

    private boolean isValidHit(Peptide p) {
        int start = p.getStart();
        int end = p.getEnd();
        if((end-start) > minPeptideLength) {
            // check enzyme specificity 
            if(protease == null || protease.checkEnzymeSpecificity(p.getParent(), start, end) >= enzymeSpecificity) { 
                return true;
            }
        }
        return false;

    }
    private boolean isValidHit(PeptideHit p) {
        int start = p.getStart();
        int end = p.getEnd();
        if((end-start) > minPeptideLength) {
            // check enzyme specificity 
            if(protease == null || protease.checkEnzymeSpecificity(p.getParent(), start, end) >= enzymeSpecificity) { 
                return true;
            }
        }
        return false;

    }
    /**

     */

    public void addFinalPeptideHit(PeptideHit p,int [] masses) {
        //numPeaks += (p.getLength()-1) * 2;
        int start = p.getStart();
        int end = p.getEnd();
        //PeptideHit p = new PeptideHit(f, start, end);
        //int numMatched = ppl.calcNumPeaksMatched(p);

        p.setProcessedPeakList(ppl);
        p.libraryCalcNumPeaksMatched(masses);

        int numMatched = p.getNumPeaksMatched();
        numPeptidesMatched++;
        //numPeaks += counts[1];
        numPeaks += p.getNumPeaks();
        totalPeptideLength += (end - start + 1);
        numPeaksMatched += p.getNumPeaksMatched();
        if (numMatched >= minNumPeaksMatched) {
            //    addPeptideHit(p);
            double prob = p.getProbability();
            int lastIndex = hits.size() - 1;
            int index = findIndex(prob, 0, lastIndex);
            if(index < NUMSCORED) {
                hits.add(index, p);
                //p.setProcessedPeakList(ppl);
                if(hits.size() > NUMSCORED) {
                    hits.remove(hits.size()-1);
                }
            }

        }
        //System.out.println("new peptide added: " + p.getSequence() + ": " + p.getNumPeaksMatched());
    }

    private void addFinalPeptideHit(PeptideHit p) {
        //numPeaks += (p.getLength()-1) * 2; 
        int start = p.getStart();
        int end = p.getEnd();
                //PeptideHit p = new PeptideHit(f, start, end);
                //int numMatched = ppl.calcNumPeaksMatched(p);

        p.setProcessedPeakList(ppl);
        p.calcNumPeaksMatched();
             
        int numMatched = p.getNumPeaksMatched();
        numPeptidesMatched++;
        //numPeaks += counts[1]; 
        numPeaks += p.getNumPeaks(); 
        totalPeptideLength += (end - start + 1);
        numPeaksMatched += p.getNumPeaksMatched();
        if (numMatched >= minNumPeaksMatched) {
        //    addPeptideHit(p);
            double prob = p.getProbability();
            int lastIndex = hits.size() - 1;
            int index = findIndex(prob, 0, lastIndex); 
            if(index < NUMSCORED) {
                hits.add(index, p);
                //p.setProcessedPeakList(ppl);
                if(hits.size() > NUMSCORED) {
                    hits.remove(hits.size()-1);
                }
            }
       
        }
        //System.out.println("new peptide added: " + p.getSequence() + ": " + p.getNumPeaksMatched()); 
    }
    public int getNumPeptideHit() {
        return hits.size();
    }
    public ArrayList<ScoredPeptideHit> getTopHits(int numTopHits) {
        /*
        if (!isSorted) {
            Collections.sort(hits);
            isSorted = true;
        }
        */
        // group PeptideHits by sequence string
        HashMap<String, ScoredPeptideHit> topHits =
                    new HashMap<String, ScoredPeptideHit>(numTopHits);
        ArrayList<ScoredPeptideHit> topHitList = new ArrayList<ScoredPeptideHit>(numTopHits);
        for (int i = 0; i < hits.size(); i++) {
            PeptideHit ph = hits.get(i);
            //String seq = ph.getExtendedSequence();
            //String seq = ph.getSequence();
            String seq = ph.getExactSequence();
//if(seq.indexOf("DDIAALVVD") != -1)
//System.out.println(seq + "\t" + hits.size());
            ScoredPeptideHit s = topHits.get(seq);
            if(s == null) {
                s = new ScoredPeptideHit(seq, 0, secondaryScoreType);
                s.setPScore(ph.getProbability());
                topHitList.add(s);
                topHits.put(seq, s);

            }
            s.addPeptideHit(ph);
            if(topHitList.size() == numTopHits) {
                break;
            }
        }
        //return topHits.values(); 
        return topHitList; 
    }
    public ArrayList<ScoredPeptideHit> getTopHits(int numTopHits, int [][] theorMassArray, boolean [] modifiedStatusArr, String[] seqArray ) {
        /*
        if (!isSorted) {
            Collections.sort(hits);
            isSorted = true;
        }
        */
        // group PeptideHits by sequence string
        HashMap<String, ScoredPeptideHit> topHits =
                new HashMap<String, ScoredPeptideHit>(numTopHits);
        ArrayList<ScoredPeptideHit> topHitList = new ArrayList<ScoredPeptideHit>(numTopHits);
        for (int i = 0; i < hits.size(); i++) {
            PeptideHit ph = hits.get(i);
            //String seq = ph.getExtendedSequence();
            //String seq = ph.getSequence();
            String seq = ph.getExactSequence();
//if(seq.indexOf("DDIAALVVD") != -1)
//System.out.println(seq + "\t" + hits.size());
            ScoredPeptideHit s = topHits.get(seq);
            if(s == null) {
                s = new ScoredPeptideHit(seq, primaryScoreType, secondaryScoreType);
                topHitList.add(s);
                topHits.put(seq, s);
            }
            s.addPeptideHit(ph);
            if(topHitList.size() == numTopHits) {
                break;
            }
        }
        for(int i=0; i<topHitList.size(); i++)
        {
            theorMassArray[i] = topHitList.get(i).getTheorMasses();
            modifiedStatusArr[i] = topHitList.get(i).isModified();
            seqArray[i] = topHitList.get(i).getSequence();
        }
        //return topHits.values();
        return topHitList;
    }

    // using binary search
    private int findIndex(double p, int low, int high) {
        int lastIndex = hits.size() - 1;
        if(lastIndex == -1 || p <= hits.get(0).getProbability()) {
            return 0;
        } else if(high == low || p >= hits.get(lastIndex).getProbability()) {
            return lastIndex + 1;
        }
        
        int mid = 0;
        while (high > low) {
            //mid = low + (high - low) / 2;
            mid = (high + low) / 2;
            double midP = hits.get(mid).getProbability();
            if (p < midP)
                high = mid;
            else if (p > midP)
                low = mid;
            else
                return mid;
            if((high-low) == 1) {
               if(p <= hits.get(low).getProbability()) {
                   return low;
               }
               if(p > hits.get(high).getProbability()) {
                   return high + 1;
               } else {
                   return high;
               }
           }
        }
        return mid;
    }

    public void setFinalResultsForPrelimScores()
    {
        scoredHits = getTopHits(20);

    }

    public void calcScoreTest()
    {
        CalcScore calcScore = new CalcScore();

        float [] xCorrResults = new float[NUMSCORED];
        int [][] theorMassArray = new int[NUMSCORED][];
        boolean [] isModifiedArray = new boolean[NUMSCORED];
        String[] seqArray = new String[NUMSCORED];
        scoredHits = getTopHits(NUMSCORED,theorMassArray,isModifiedArray,seqArray);


        int [] xcorrRank = new int[NUMSCORED];
        float [] massCorr = ppl.getMassCorrFloat();
        float [] massSum = ppl.getMassSumFloat();
        float [] result = new float[2];
        printTheoreticalMass(theorMassArray);
        printMassCorr(massCorr,massSum);
        calcScore.calcScore(xCorrResults,xcorrRank,result,theorMassArray,massCorr,massSum, ppl.ACCURACYFACTOR,seqArray, isModifiedArray);
        printXCorrResults(xCorrResults);
        primaryScoreMean = result[0];
        primaryScoreDeviation = result[1];
        int i=0;
        for(ScoredPeptideHit p : scoredHits) {
            //System.out.println(p.getXCorr());
            p.setXCorr(xCorrResults[i]);
            p.setXCorrRank(xcorrRank[i]);
            double score = -Math.log10(DistributionCalculator.getBinomialSum(ppl.getPTrue(), p.getNumPeaks(), p.getNumPeaksMatched()));

            p.setPScore(score);
            double zscore = (p.getXCorr() - primaryScoreMean)/primaryScoreDeviation;
            double pvalue = StatCalc.zScore2PValue(zscore);
            double evalue = pvalue*numPeptidesMatched;
            p.setZscore(zscore);

            p.setPValue(pvalue);
            p.setEValue(evalue);
            i++;
//System.out.println("Z\t" + zscore + "\tpvalue\t" + pvalue + "\tevalue\t" + evalue);
            //p.setZscore(((p.getXCorr() - primaryScoreMean)/primaryScoreDeviation)/2); //temperary change for CAMSI
            //p.setXCorr((p.getPrimaryScore() - primaryScoreMean)/primaryScoreDeviation);
        }
        Collections.sort(scoredHits);
        finalResult = getFinalResult();
    }
    public void printMassCorr(float [] massCorr, float [] massSum)
    {
        System.out.println("Mass Corr");
        System.out.println(Arrays.toString(massCorr));
        System.out.println("Mass Sum");
        System.out.println(Arrays.toString(massSum));

    }

    public void printTheoreticalMass(int [][] theorMass)
    {
        System.out.println("Theoretical Mass");
        for(int i=0; i<theorMass.length; i++)
        {
            System.out.println(Arrays.toString(theorMass[i]));
        }
    }
    public void printXCorrResults(float [] xcorrResults)
    {
        System.out.println("Xcorr results");
        System.out.println(Arrays.toString(xcorrResults));


    }

    public void addScoredPeptideHit(ScoredPeptideHit sph)
    {

    }

    public void setPrelimScoreHits(List<ScoredPeptideHit> sphList, int size)
    {
        Collections.sort(sphList, new Comparator<ScoredPeptideHit>() {
            @Override
            public int compare(ScoredPeptideHit scoredPeptideHit, ScoredPeptideHit t1) {
                return Double.compare(scoredPeptideHit.getPScore(),t1.getPScore());
            }
        });
        finalResult = sphList.subList(0,size);
        numPeptidesMatched = finalResult.size();
    }

    public static int getNUMFINALRESULT()
    {
        return NUMFINALRESULT;
    }

    public String getMs2Name() {
        return ms2Name;
    }

    public void setMs2Name(String ms2Name) {
        this.ms2Name = ms2Name;
    }

    public double getRetTime() {
        return retTime;
    }

    public void setRetTime(double retTime) {
        this.retTime = retTime;
    }
}
