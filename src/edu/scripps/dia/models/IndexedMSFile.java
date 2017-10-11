package edu.scripps.dia.models;

import edu.scripps.dia.util.MSIndexBuilder;
import gnu.trove.TFloatArrayList;
import gnu.trove.TIntArrayList;
import gnu.trove.TIntLongHashMap;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Map;

/**
 * Created by yateslab on 7/20/17.
 */
public class IndexedMSFile implements AutoCloseable{

    private TIntLongHashMap scanPosMap;
    private String indexedFilePath;
    private String msFilePath;
    private RandomAccessFile file;

    public IndexedMSFile(String msFile) throws IOException
    {

        msFilePath = msFile+".ms2";
        indexedFilePath = msFilePath+".index";
        init();
    }



    public IndexedMSFile(String indexedFile, String msFile) throws IOException
    {
        indexedFilePath = indexedFile;
        msFilePath = msFile;
        init();
    }

    private void init() throws IOException {
        MSIndexBuilder builder = new MSIndexBuilder(indexedFilePath);
        scanPosMap = builder.readMS2IndexFile();
        file = new RandomAccessFile(msFilePath,"r");
    }

    public void getSpectra(int scan, int[] consensusSpectra) throws IOException {
        goToScan(scan);
        String line;
        while((line = readLine()) !=null)
        {
            if(Character.isDigit(line.charAt(0)))
            {
                String [] arr = line.split("\\s");
            }
        }

    }


    @Override
    public void close() throws Exception {
        file.close();
    }

    public void goToScan(int scan) throws IOException {
        long pos = scanPosMap.get(scan);
        file.seek(pos);
    }

    public String readLine() throws IOException {
       return file.readLine();
    }

    public String readLine(int scan) throws IOException {
        goToScan(scan);
        return file.readLine();
    }

    public TIntLongHashMap getScanPosMap() {
        return scanPosMap;
    }

    public String getIndexedFilePath() {
        return indexedFilePath;
    }

    public String getMsFilePath() {
        return msFilePath;
    }

    public RandomAccessFile getFile() {
        return file;
    }
}
