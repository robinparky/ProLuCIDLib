package edu.scripps.dia;

import edu.scripps.pms.mspid.ScoredPeptideHit;

import java.util.Comparator;

/**
 * Created by yateslab on 10/12/17.
 */
public class ScoredPeptideHitComparator implements Comparator<ScoredPeptideHit> {
    @Override
    public int compare(ScoredPeptideHit scoredPeptideHit, ScoredPeptideHit t1) {
        return Double.compare(scoredPeptideHit.getPScore(),t1.getPScore());
    }


}
