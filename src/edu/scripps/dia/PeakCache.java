package edu.scripps.dia;

import edu.scripps.pms.mspid.ProcessedPeakList;

import org.apache.commons.collections4.map.LRUMap;


/**
 * Created by Titus Jung titusj@scripps.edu on 12/14/17.
 */
public class PeakCache extends LRUMap<Integer, ProcessedPeakList> {

    private static final int MAX_ENTRIES = 500;
    public PeakCache()
    {
        super(MAX_ENTRIES+1,0.75f,true);
    }

    protected boolean removeLRU(LinkEntry<Integer, ProcessedPeakList> entry)
    {
        entry.getValue().dumpBoolMass();
        return true;
    }

}
