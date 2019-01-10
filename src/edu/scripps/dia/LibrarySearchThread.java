package edu.scripps.dia;

import edu.scripps.pms.mspid.SearchResult;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;

import java.util.Iterator;

/**
 * Created by Titus Jung titusj@scripps.edu on 4/11/18.
 */
public class LibrarySearchThread implements Runnable{

    private LibrarySearchEngine lse;
    private String msName;

    public LibrarySearchThread(LibrarySearchEngine lse, String msName) {
        this.lse = lse;
        this.msName = msName;
    }

    @Override
    public void run() {
        PeakList sourceList;
        try
        {
            while((sourceList=lse.getPeakList())!=null)
            {
                for (Iterator<Zline> it = sourceList.getZlines(); it.hasNext(); ) {
                    Zline zline = it.next();
                    SearchResult r = lse.searchIndexedDB(sourceList,zline);
                    r.setMs2Name(msName);
                    lse.write(r.outputResults());
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
