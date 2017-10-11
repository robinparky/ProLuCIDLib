package edu.scripps.dia.util;

import gnu.trove.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version $Id: MSIndexBuilder.java,v 1.6 2014/01/31 21:47:14 rpark Exp $
 */
public class MSIndexBuilder
{
    private TIntLongHashMap posMap = new TIntLongHashMap();
    private Hashtable<Integer, String> indexHt = new Hashtable<Integer, String>();
    private TIntDoubleHashMap precursorMap = new TIntDoubleHashMap();
    private TDoubleDoubleHashMap rtPrecursorMap = new TDoubleDoubleHashMap();
    private TDoubleLongHashMap rtPosMap = new TDoubleLongHashMap();
    private TDoubleIntHashMap rtScanNumMap = new TDoubleIntHashMap();
    private Hashtable<Integer, String> scanTypeHt = new Hashtable<Integer, String>();
    private Hashtable<Integer, Integer> precScanMap = new Hashtable<Integer, Integer>();
    
    //private RandomAccessFile file;
    private File indexedFile;
    private String fileName;

    public static void main(String args[]) throws Exception
    {
            //MSIndexBuilder builder = new MSIndexBuilder("/hdb/relex/data/RelEx_new_data/data-independent/100803_yeast_1to1_DI_FIM_LTQ-01.ms1");
            //builder.buildIndex();

    }

    public MSIndexBuilder(String fileName)
    {
        this.indexedFile = new File(fileName);
    }

    public MSIndexBuilder(File indexedFile, String fileName)
    {
        this.indexedFile = indexedFile;
        this.fileName = fileName;

    }
/*
    public MSIndexBuilder(RandomAccessFile file)
    {
        this.file = file;
    }
    
    public MSIndexBuilder(File indexedFile, RandomAccessFile file)
    {
        this.file = file;
        this.indexedFile = indexedFile;
    }
*/
    public TIntDoubleHashMap getPrecursorMap()
    {
        return precursorMap;
    }

    public TIntLongHashMap readIndexFile() throws IOException
    {

        if(fileName.endsWith("ms1") || fileName.endsWith("mszm"))
            return readMS1IndexFile();
        else if( fileName.endsWith("ms2") )
            return readMS2SimpleIndexFile();
	else if(fileName.endsWith("ms2") || fileName.endsWith("ms3"))
            return readMS2IndexFile();
        else
            throw new IOException("File extension is incorrect");
           
    }

    public TIntLongHashMap readMS1IndexFile() throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(indexedFile));
        String lastLine;

        while( (lastLine=br.readLine()) != null)
        {
            String[] str = lastLine.split("\t");
            posMap.put(Integer.parseInt(str[0]), Long.parseLong(str[1]));

        //    System.out.println(" ===>>" + Integer.parseInt(str[0]) + " " + Long.parseLong(str[1]));
        }
    
        return posMap;
    }

    public TIntLongHashMap readMS2SimpleIndexFile() throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(indexedFile));
        String curLine;

        String lastLine;
        
        while( (lastLine=br.readLine()) != null)
        {
            String[] str = lastLine.split("\t");
            
            int sNum = Integer.parseInt(str[0]);
            long pos = Long.parseLong(str[1]);

	    posMap.put(sNum, pos);
        }
    
        return posMap;

    }

    public TIntLongHashMap readMS2IndexFile() throws IOException
    {
        
        String lastLine;
        BufferedReader br = new BufferedReader(new FileReader(indexedFile));
        
        while( (lastLine=br.readLine()) != null)
        {
            String[] str = lastLine.split("\t");
            int scanNum = Integer.parseInt(str[0]);
            long pos = Long.parseLong(str[1]);
            double pCursor = Double.parseDouble(str[2]);
            double rt;

	    if(str.length<4)
		rt = -1.0;
	    else
		rt = Double.parseDouble(str[3]);
                        
            posMap.put( scanNum, pos);
	    if(str.length>=4) {
		indexHt.put(scanNum, lastLine);
                scanTypeHt.put(scanNum, str[4]);
            }

            precursorMap.put( scanNum, pCursor );            
            rtPrecursorMap.put(rt, pCursor);
            rtScanNumMap.put(rt, scanNum);
            
	    if(str.length>3)
		rtPosMap.put( rt, pos );
            
            
            if(str.length>5)
                precScanMap.put(Integer.parseInt(str[5]), scanNum);
                
            
            
            
        }
        
        
        
        return posMap;
    }

    public TDoubleLongHashMap getRtPosMap() {
        return rtPosMap;
    }

    public void setRtPosMap(TDoubleLongHashMap rtPosMap) {
        this.rtPosMap = rtPosMap;
    }

    public TDoubleDoubleHashMap getRtPrecursorMap() {
        return rtPrecursorMap;
    }

    public void setRtPrecursorMap(TDoubleDoubleHashMap rtPrecursorMap) {
        this.rtPrecursorMap = rtPrecursorMap;
    }

    public TDoubleIntHashMap getRtScanNumMap() {
        return rtScanNumMap;
    }

    public void setRtScanNumMap(TDoubleIntHashMap rtScanNumMap) {
        this.rtScanNumMap = rtScanNumMap;
    }

    public Hashtable<Integer, String> getIndexHt() {
        return indexHt;
    }

    public void setIndexHt(Hashtable<Integer, String> indexHt) {
        this.indexHt = indexHt;
    }
    
    public String getScanType(int scanNum) {
        return this.scanTypeHt.get(scanNum);
    }

    public Hashtable<Integer, String> getScanTypeHt() {
        return scanTypeHt;
    }

    public void setScanTypeHt(Hashtable<Integer, String> scanTypeHt) {
        this.scanTypeHt = scanTypeHt;
    }

    public Hashtable<Integer, Integer> getPrecScanMap() {
        return precScanMap;
    }

    public void setPrecScanMap(Hashtable<Integer, Integer> precScanMap) {
        this.precScanMap = precScanMap;
    }

    
            
}
