package edu.scripps.dia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Titus Jung titusj@scripps.edu on 5/15/18.
 */
public class LibraryCSV {

    public static long peptideID = 0;
    public static void main(String [] args) throws  IOException {
        String libAddress = args[0];
        File f = new File(libAddress);
        if(f.exists()) f.delete();
        LibraryIndexer libraryIndexer = new LibraryIndexer();

        Map<String,Long> seqKeyMap = new HashMap<>();
        for(int i=1; i<args.length; i++)
        {
            System.out.println("reading "+args[i]);
            String path = args[i].substring(0,args[i].lastIndexOf(File.separatorChar));

            libraryIndexer.readDTASelectCreateCSV(args[i],path,libAddress,null,seqKeyMap);
        }
    }

}
