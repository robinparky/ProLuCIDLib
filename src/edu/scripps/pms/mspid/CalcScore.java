package edu.scripps.pms.mspid;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by rpark on 5/1/17.
 */
public class CalcScore {
    public class PeptideHitStruct
    {
        int id;
        float xcorr;
        int rank;
        String seq;
        boolean isModified;
    }

    class PeptideHitXCorrComparator implements Comparator<PeptideHitStruct>
    {
        @Override
        public int compare(PeptideHitStruct o1, PeptideHitStruct o2) {
            int result=  Float.compare(o1.xcorr,o2.xcorr);
            if(result!=0) return -result;
            else
            {
                if(o1.isModified==o2.isModified) return 0;
                else if(o2.isModified) return -1;
                else return 1;
            }
        }
    }
    class PeptideHitIDComparator implements Comparator<PeptideHitStruct>
    {
        @Override
        public int compare(PeptideHitStruct o1, PeptideHitStruct o2) {
            return o1.id-o2.id;
        }
    }


    public double calcXCorr(int[] theorMass, double [] massCorr,double [] massSum, double ACCURACYFACTOR ) {
        int maxIndex = massCorr.length - 75;
        double sumProduct = 0;
        double newMoreSumProduct = 0;
        for(int i=0; i<theorMass.length;i++) {
            if (theorMass[i] > 0 && i > 75 && i < maxIndex) {
                int mass = i;
                double intensity = (int) (theorMass[i] * ACCURACYFACTOR);
                sumProduct += massCorr[mass] * intensity;
                newMoreSumProduct += massSum[mass] * intensity;
            }
        }
        double xcorr = (0.993377483f * sumProduct - newMoreSumProduct / 151) / 10000;
        xcorr = xcorr < 0.00001 ? 0.00001f : xcorr;
        return xcorr;
    }
    public  void calcScore(float [] xscorrResults, int[] xcorrRank, float [] resultArr, int[][] theorMass, float [] massCorr,
                           float [] massSum, float ACCURACYFACTOR, String[] seqArray, boolean [] modifiedArray )
    {
        int maxIndex = massCorr.length - 75;
        for(int i=0; i<theorMass.length; i++)
        {
            float sumProduct = 0;
            float newMoreSumProduct = 0;
            for(int j=0; j<theorMass[i].length;j++) {
                if (theorMass[i][j] > 0 && j > 75 && j < maxIndex) {
                    int mass = j;
                    double intensity = (int) (theorMass[i][j] * ACCURACYFACTOR);
                    sumProduct += massCorr[mass] * intensity;
                    newMoreSumProduct += massSum[mass] * intensity;
                }
            }
            float xcorr = (0.993377483f * sumProduct - newMoreSumProduct / 151) / 10000;
            xcorr = xcorr < 0.00001 ? 0.00001f : xcorr;
            xscorrResults[i] = xcorr;
        }
        float [] xcorrCopy = Arrays.copyOf(xscorrResults,xscorrResults.length);
        PeptideHitStruct [] phArray = new PeptideHitStruct[xscorrResults.length];
        for(int i=0; i<xscorrResults.length; i++)
        {
            PeptideHitStruct phs = new PeptideHitStruct();
            phs.id =i;
            phs.xcorr = xscorrResults[i];
            phs.seq = seqArray[i];
            phs.isModified = modifiedArray[i];
            phArray[i] = phs;
        }
        Arrays.sort(phArray,new PeptideHitXCorrComparator());
        int rank = 0;
        double lastScore = 100;
        for(int i=0; i<phArray.length; i++)
        {
            if(phArray[i].xcorr<lastScore)
            {
                rank++;
                lastScore = phArray[i].xcorr;
            }
            phArray[i].rank = rank;
        }
        String topseq = null;
        float primaryScoreMean = 0;
        int numElements =0;
        for(int i=0; i<phArray.length; i++)
        {
            if(topseq == null)
            {
                topseq = phArray[i].seq;
            }
            if(phArray[i].rank != 1 && (!phArray[i].seq.equals(topseq)))
            {
                primaryScoreMean += phArray[i].xcorr;
                numElements++;
            }
        }
        float primaryScoreDeviation =0;
        if(numElements > 2)
        {
            primaryScoreMean /=numElements;
            for(int i=0; i<phArray.length; i++)
            {

                if(phArray[i].rank != 1 && (!phArray[i].seq.equals(topseq)))
                {
                    float diff = phArray[i].xcorr - primaryScoreMean;
                    primaryScoreDeviation += diff *diff;
                }
            }
            primaryScoreDeviation =(float) Math.sqrt((primaryScoreDeviation/(numElements-1)));
            if(primaryScoreDeviation ==0)
            {
                primaryScoreDeviation = 1;
            }
        }
        Arrays.sort(phArray,new PeptideHitIDComparator());
        for(int i=0; i<phArray.length; i++)
        {
            xcorrRank[i] = phArray[i].rank;

        }

        resultArr[0] = primaryScoreMean;
        resultArr[1] = primaryScoreDeviation;
    }



}
