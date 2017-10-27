package edu.scripps.dia.models;

import gnu.trove.TIntLongHashMap;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Map;

/**
 * Created by Titus Jung titusj@scripps.edu on 10/27/17.
 */
public class ZIndexedMSFile implements AutoCloseable {

    private Map<Integer,TIntLongHashMap> zScanPosMap;
    private String indexedFilePath;
    private String msFilePath;
    private RandomAccessFile file;

    public ZIndexedMSFile(String msFile) throws IOException
    {

        msFilePath = msFile+".ms2";
        indexedFilePath = msFilePath+".zindex";
    }

    private void init()
    {

    }


    @Override
    public void close() throws Exception {
        if(file!=null)
        {
            file.close();
        }
    }
}
