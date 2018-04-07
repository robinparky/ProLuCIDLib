/*
 * IndexedFile.java
 *
 * Created on March 21, 2005, 11:51 AM
 */

package edu.scripps.dia.util;



//import edu.scripps.pms.censu.util.io.MzxmlSpectrumReader;
import gnu.trove.*;

import java.io.*;
import java.util.*;


/**
 *
 * @author rpark
 */
public class IndexedFile {

    /** Creates a new instance of IndexedFile */
    private String fileName;

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    private RandomAccessFile file=null;
    private MzxmlSpectrumReader mzreader = null;
    private int[] keys;
    private double[] rtArr;
    private TIntLongHashMap index;
    private MSIndexBuilder builder;
    private File indexFile;
    private int startScanNum;
    private int lastScanNum;
    private TIntDoubleHashMap precursorMap;
    private TIntDoubleHashMap retentionTimeMap;//scan number to retention time
    private TIntDoubleHashMap ionInjectionMap;//scan number to ionInjection time
    private TDoubleIntHashMap retentonToScanMap = new TDoubleIntHashMap(); // retention Time to scan NUmebr
    private double[] retentionTimeSortedArr = null;
    private TIntArrayList keyIntArr;
    private TDoubleLongHashMap rtPosMap;
    private TIntLongHashMap scanPositionMap;
    private TDoubleIntHashMap rtScanNumMap;
    private TDoubleDoubleHashMap rtPrecursorMap = new TDoubleDoubleHashMap();
    private Hashtable<Integer, String> indexHt = new Hashtable<Integer, String>();
    private Hashtable<Integer, String> scanTypeHt = new Hashtable<Integer, String>();
    private TIntIntHashMap idScanToMs1ScanMap = new TIntIntHashMap();

    private Hashtable<Integer, Integer> precScanMap = new Hashtable<Integer, Integer>();
    private TIntArrayList scanList = new TIntArrayList();
    //    private int numOfIsoWindows=-1;
//    private double[] windowStartValues;
    //private double isolationWindow;
    private final int FULL_SCAN = 1;
    private TIntDoubleHashMap scanRtMap = new TIntDoubleHashMap();
    private TIntArrayList scanArr = new TIntArrayList();
    private TIntIntHashMap scanToIndexMap = new TIntIntHashMap();

    public static void main(String args[]) throws Exception
    {

        // if(true) return;

        File f = new File(args[0]);

        IndexedFile file = new IndexedFile(f, args[1]);
    }

    //mzxml file
    public IndexedFile(String fileName) throws Exception
    {
        this.fileName = fileName;

        createMzXMLScans();


    }


    public void createMzXMLScans() throws Exception
    {
        keyIntArr = new TIntArrayList();
        this.mzreader = new MzxmlSpectrumReader(fileName);
        for(Iterator<String> itr=mzreader.getScanNums(FULL_SCAN); itr.hasNext(); )
        {
            String scanNum = itr.next();
            keyIntArr.add( Integer.parseInt(scanNum) );
        }

        keys = keyIntArr.toNativeArray();
    }

    public IndexedFile(File indexFile, String fileName) throws IOException
    {
        this.indexFile = indexFile;
        if(!indexFile.exists()) {
            MSIndexFileCreator.createIndexFile(fileName);
        }

        this.fileName = fileName;
        readIndexFile();
        createScanNum();

    }

    public IndexedFile(String indexFile, String fileName) throws IOException
    {
        this.indexFile = new File(indexFile);
        if(!this.indexFile.exists()) {
            MSIndexFileCreator.createIndexFile(fileName);
        }

        this.fileName = fileName;
        readIndexFile();
        createScanNum();

    }

    private void createScanNum()
    {
        startScanNum = keys[0];
        lastScanNum = keys[keys.length-1];
    }

    private void readIndexFile() throws IOException
    {
        file = new RandomAccessFile(fileName, "r");

        builder = new MSIndexBuilder(indexFile, fileName);
        index = builder.readIndexFile();
        retentionTimeMap = builder.getRetentionMap();
        setIonInjectionMap(builder.getIonInjectionMap());
        retentonToScanMap = builder.getRetToScanNumerMap();
        precursorMap = builder.getPrecursorMap();
        this.rtPrecursorMap = builder.getRtPrecursorMap();
        rtPosMap = builder.getRtPosMap();
        scanPositionMap = builder.getPosMap();
        rtScanNumMap = builder.getRtScanNumMap();
        keys = index.keys();
        //rtArr = rtPosMap.keys();//robin check this value again
        rtArr = builder.getRetArr().toNativeArray();
        indexHt = builder.getIndexHt();
        scanArr = builder.getScanArr();
        scanTypeHt = builder.getScanTypeHt();
        precScanMap = builder.getPrecScanMap();
        idScanToMs1ScanMap = builder.getIdScanToMs1ScanMap();
        scanList = builder.getScanList();
        this.scanRtMap = builder.getScanRtMap();
        this.scanToIndexMap = builder.getScanToIndexMap();

        Arrays.sort(keys);
        file.close();
        //  Arrays.sort(rtArr);
    }

   /*
    private void createIndexFile() throws IOException
    {
        file = new RandomAccessFile(fileName, "r");

        builder = new MSIndexBuilder(file);
        index = builder.buildIndex();

        keys = index.keys();
        Arrays.sort(keys);
    }
    */

    public void close() throws IOException
    {
        if(null != file)
            file.close();
    }

    public TIntLongHashMap getMsIndex()
    {
        return index;
    }

    public TIntDoubleHashMap getPrecursorMap()
    {
        return precursorMap;
    }

    /**
     * returns the Scan Number array
     * @return
     */
    public int[] getKeys()
    {
        return keys;
    }

    public int getStartScanNum()
    {
        return startScanNum;
    }

    public int getLastScanNum()
    {
        return lastScanNum;
    }

    public RandomAccessFile getFile()
    {
        //check if file is closed or not
        try {
            file.length();
        } catch(Exception e) {
            try {
                file = new RandomAccessFile(fileName, "r");
            } catch (IOException ioe) {}
        }

        return file;
    }

   /*
    public int getNumOfIsoWindows()
    {
        return numOfIsoWindows;
    }
/*
    public double getIsolationWindow() {
        return isolationWindow;
    }
*/

    public MzxmlSpectrumReader getMzreader() {
        return mzreader;
    }

    public void setMzreader(MzxmlSpectrumReader mzreader) {
        this.mzreader = mzreader;
    }

    public TIntArrayList getKeyIntArr() {
        return keyIntArr;
    }

    public void setKeyIntArr(TIntArrayList keyIntArr) {
        this.keyIntArr = keyIntArr;
    }

    public long getPosByRt(double d) {
        return this.rtPosMap.get(d);
    }

    public TDoubleLongHashMap getRtPosMap() {
        return rtPosMap;
    }

    public void setRtPosMap(TDoubleLongHashMap rtPosMap) {
        this.rtPosMap = rtPosMap;
    }

    public double[] getRtArr() {
        return rtArr;
    }

    public void setRtArr(double[] rtArr) {
        this.rtArr = rtArr;
    }

    public TDoubleDoubleHashMap getRtPrecursorMap() {
        return rtPrecursorMap;
    }

    public void setRtPrecursorMap(TDoubleDoubleHashMap rtPrecursorMap) {
        this.rtPrecursorMap = rtPrecursorMap;
    }

    public int getScanNumByRtIndex(int rtIndex)
    {
        return this.rtScanNumMap.get( this.rtArr[rtIndex] );
    }

    public TDoubleIntHashMap getRtScanNumMap() {
        return rtScanNumMap;
    }

    public void setRtScanNumMap(TDoubleIntHashMap rtScanNumMap) {
        this.rtScanNumMap = rtScanNumMap;
    }

    public double getPrecursorByScanNum(int scanNum)
    {

        return this.precursorMap.get(scanNum);
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

    public int getPrecursorScan(int daughterScan) {

        Integer into = this.precScanMap.get(daughterScan);
        if(null == into) return -1;

        return into;
        //return this.precScanMap.get(daughterScan);
    }

    /**
     * ScanNumber to Retention Time Map
     * @return the retentionTimeMap
     */
    public TIntDoubleHashMap getRetentionTimeMap() {
        return retentionTimeMap;
    }

    /**
     * @param retentionTimeMap the retentionTimeMap to set
     */
    public void setRetentionTimeMap(TIntDoubleHashMap retentionTimeMap) {
        this.retentionTimeMap = retentionTimeMap;
    }

    /**
     * @return the ionInjectionMap
     */
    public TIntDoubleHashMap getIonInjectionMap() {
        return ionInjectionMap;
    }

    /**
     * @param ionInjectionMap the ionInjectionMap to set
     */
    public void setIonInjectionMap(TIntDoubleHashMap ionInjectionMap) {
        this.ionInjectionMap = ionInjectionMap;
    }

    /**
     *
     * @return Retentiontime to ScanNumber hash map....
     */
    public TDoubleIntHashMap getRetentonToScanMap() {
        return retentonToScanMap;
    }

    public void setRetentonToScanMap(TDoubleIntHashMap retentonToScanMap) {
        this.retentonToScanMap = retentonToScanMap;
    }

    public TIntArrayList getScanList() {
        return scanList;
    }

    public void setScanList(TIntArrayList scanList) {
        this.scanList = scanList;
    }

    public TIntDoubleHashMap getScanRtMap() {
        return scanRtMap;
    }

    public void setScanRtMap(TIntDoubleHashMap scanRtMap) {
        this.scanRtMap = scanRtMap;
    }



    /**
     * returns sorted key of retentionTime..
     * @return
     */
    public double[] getRetentionTimeSortedArr() {
        if(retentionTimeSortedArr==null)
        {
            retentionTimeSortedArr = retentonToScanMap.keys();
            Arrays.sort(retentionTimeSortedArr);
        }
        return retentionTimeSortedArr;
    }

    public void setRetentionTimeSortedArr(double[] retentionTimeSortedArr) {
        this.retentionTimeSortedArr = retentionTimeSortedArr;
    }

    /**
     * Gives the index of the scanNumber from the keys array
     * @param scanNumber
     */
    //don't use this.  this is very expensive
    public int getIndexFromScanNumber(int scanNumber)
    {
        System.out.println("IndexedFile.getIndexFormScanNumber: warning.. this method is expensive");
        int index = Arrays.binarySearch(keys,scanNumber);
        if(index<0) //Cannot find index
            index=-(++index); //Math.abs(++keyIndex);

        if(index>=keys.length)
            index--;
        return index;
    }

    public TIntLongHashMap getScanPositionMap() {
        return scanPositionMap;
    }

    public void setScanPositionMap(TIntLongHashMap scanPositionMap) {
        this.scanPositionMap = scanPositionMap;
    }

    public long getPositionByScan(int scan) {
        return this.scanPositionMap.get(scan);
    }

    public long getPositionByIndex(int index) {
        return getPositionByScan( this.scanArr.get(index) );

    }

    public TIntArrayList getScanArr() {
        return scanArr;
    }

    public void setScanArr(TIntArrayList scanArr) {
        this.scanArr = scanArr;
    }

    public TIntIntHashMap getScanToIndexMap() {
        return scanToIndexMap;
    }

    public void setScanToIndexMap(TIntIntHashMap scanToIndexMap) {
        this.scanToIndexMap = scanToIndexMap;
    }

    public int getIndexByScan(int scanNum) {
        return this.scanToIndexMap.get(scanNum);
    }

    public TIntIntHashMap getIdScanToMs1ScanMap() {
        return idScanToMs1ScanMap;
    }

    public void setIdScanToMs1ScanMap(TIntIntHashMap idScanToMs1ScanMap) {
        this.idScanToMs1ScanMap = idScanToMs1ScanMap;
    }

    public File getIndexFile() {
        return indexFile;
    }
}
