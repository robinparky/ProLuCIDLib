package edu.scripps.dia;

import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.mspid.Modifications;
import edu.scripps.pms.mspid.ProcessedPeakList;
import gnu.trove.TIntArrayList;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by yateslab on 8/17/17.
 */
public class MassRangeFinder {

    public static boolean withinRange(double source, double target, SearchParams params)
    {
        if(helper(source,target,params)) return true;
        for(Iterator<Modifications> it = params.getAllModifications(); it.hasNext();) {
            Modifications m = it.next();
            if(m != null && m.getDiffModsShift() != 0 ) {
                if(helper(source-(float)m.getMassShift(),target, params)) return true;
            }
        }
        return false;
    }


    private static boolean helper(double source, double target, SearchParams params)
    {

        double acc = params.getPrecursorTolerance()/1000000.0f;

        if(params.getNumIsotopicPeaks() == 0) { // for low resolution, traditional sequest like
            double high = ((source + params.getHighPrecursorTolerance()/1000));
            double low  = ((source - params.getLowPrecursorTolerance()/1000));
            return target<high && target>low;
        } else if(params.getNumIsotopicPeaks() == 1) { // for deisotoped high resolution data
            double diffs = source*acc;
            double high = ((source + diffs));
            double low  = ((source - diffs));
            return target<high && target>low;

        } else { //for non deisotoped high resolution data
            double diffs = source*acc*2;
            for(int i = 0; i < params.getNumIsotopicPeaks(); i++) {
                double low = source - diffs/2 - i* MassSpecConstants.MASSDIFFC12C13;
                double high = source + diffs/2 - i* MassSpecConstants.MASSDIFFC12C13;
                if(target<high && target>low) return true;
            }
        }
        return false;
    }

    public static void findRange(double prcMass, SearchParams params, TIntArrayList lowLimits, TIntArrayList highlimits)
    {
        helper(prcMass,params,lowLimits,highlimits);
        for(Iterator<Modifications> it = params.getAllModifications(); it.hasNext();) {
            Modifications m = it.next();
            if(m != null && m.getDiffModsShift() != 0 ) {
                helper(prcMass-(float)m.getMassShift(),params,lowLimits,highlimits);
            }
        }
    }

    private static void helper(double mass, SearchParams params, TIntArrayList lowLimits, TIntArrayList highLimits)
    {

        double acc = params.getPrecursorTolerance()/1000000.0f;

        if(params.getNumIsotopicPeaks() == 0) { // for low resolution, traditional sequest like
            int highLimit =(int) ((mass + params.getHighPrecursorTolerance()/1000)*1000);
            int lowLimit =(int) ((mass - params.getLowPrecursorTolerance()/1000)*1000);
            lowLimits.add(lowLimit);
            highLimits.add(highLimit);
        } else if(params.getNumIsotopicPeaks() == 1) { // for deisotoped high resolution data
            double diffs = mass*acc;
            int highLimit = (int)((mass + diffs)*1000);
            int lowLimit = (int)((mass - diffs)*1000);
            lowLimits.add(lowLimit);
            highLimits.add(highLimit);
        } else { //for non deisotoped high resolution data

            double diffs = mass*acc*2;
            for(int i = 0; i < params.getNumIsotopicPeaks(); i++) {
                lowLimits.add((int)((mass - diffs/2 - i* MassSpecConstants.MASSDIFFC12C13)*1000));
                highLimits.add((int)((lowLimits.get(i)+diffs*1000)));
            }

        }
    }

}
