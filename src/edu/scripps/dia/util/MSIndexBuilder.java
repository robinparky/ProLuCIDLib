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

    private TIntDoubleHashMap retMap = new TIntDoubleHashMap();//Scan NUmbe to Retention time
    private TIntDoubleHashMap ionInjectionMap = new TIntDoubleHashMap();// scan NUmebr to ion injection time
    private TDoubleIntHashMap retToScanNumerMap = new TDoubleIntHashMap();// retentionTime to scan NUmber
    private Hashtable<Integer, String> indexHt = new Hashtable<Integer, String>();
    private TIntDoubleHashMap precursorMap = new TIntDoubleHashMap();
    private TDoubleDoubleHashMap rtPrecursorMap = new TDoubleDoubleHashMap();
    private TDoubleLongHashMap rtPosMap = new TDoubleLongHashMap();
    private TDoubleIntHashMap rtScanNumMap = new TDoubleIntHashMap();
    private Hashtable<Integer, String> scanTypeHt = new Hashtable<Integer, String>();
    private Hashtable<Integer, Integer> precScanMap = new Hashtable<Integer, Integer>();
    private TDoubleArrayList retList = new TDoubleArrayList();
    private TIntArrayList scanList = new TIntArrayList();
    private TIntDoubleHashMap scanRtMap = new TIntDoubleHashMap();
    private TDoubleArrayList retArr = new TDoubleArrayList();
    private TIntArrayList scanArr = new TIntArrayList();
    private TIntIntHashMap scanToIndexMap = new TIntIntHashMap();
    private TIntIntHashMap idScanToMs1ScanMap = new TIntIntHashMap();

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
        int count=0;

        while( (lastLine=br.readLine()) != null)
        {
            String[] str = lastLine.split("\t");
            try {
                Integer.parseInt((str[0]));
            } catch(Exception e) {
                System.out.println("error");
            }
            int scanNum =Integer.parseInt(str[0]);
            posMap.put(scanNum, Long.parseLong(str[1]));
            this.scanToIndexMap.put(scanNum, count++);
            scanArr.add(scanNum);
            if(str.length>2)
            {
                int tmpScan = scanNum;
                double tmpRt = Double.parseDouble(str[2]);


                retMap.put(scanNum, tmpRt);
                retArr.add(tmpRt);

                ionInjectionMap.put(tmpScan, Double.parseDouble(str[3]));
                retToScanNumerMap.put(tmpRt, tmpScan);
                scanRtMap.put(tmpScan, tmpRt);
            }
            //    System.out.println(" ===>>" + Integer.parseInt(str[0]) + " " + Long.parseLong(str[1]));
        }

        return posMap;
    }

    public TDoubleArrayList getRetList() {
        return retList;
    }

    public void setRetList(TDoubleArrayList retList) {
        this.retList = retList;
    }


    public TIntDoubleHashMap getRetentionMap ()
    {
        return retMap;
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
        int count=0;

        while( (lastLine=br.readLine()) != null)
        {
            String[] str = lastLine.split("\t");
            int scanNum = Integer.parseInt(str[0]);
            this.scanToIndexMap.put(scanNum, count++);

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
            {
                precScanMap.put(Integer.parseInt(str[5]), scanNum);
                idScanToMs1ScanMap.put(scanNum,Integer.parseInt(str[5]));
            }




        }



        return posMap;
    }


/*
    public TIntLongHashMap readMS2IndexFile() throws IOException
    {
        Configuration conf = Configuration.getInstance();
        int numIsoWindows = conf.getNumOfIsolationWindow();
        double isoWin = conf.getIsolationWindow();
        double[] arr = conf.getPrecursorArr();
        double lastPrecur=-1;

        if(null != arr)
            lastPrecur = arr[arr.length-1];

        double massDiff = isoWin - conf.getMassRange(); //this is negative value

        BufferedReader br = new BufferedReader(new FileReader(indexedFile));
        String curLine;

        String lastLine = br.readLine();
        String[] str;
        str = lastLine.split("\t");

        int prevSNum = Integer.parseInt(str[0]);
        double prevPrecur = Double.parseDouble(str[2]);
        double curPrecur;
        posMap.put(prevSNum, Long.parseLong(str[1]));
        precursorMap.put(prevSNum, prevPrecur);

        while( (lastLine=br.readLine()) != null)
        {
            str = lastLine.split("\t");

            curPrecur = Double.parseDouble(str[2]);

            double curDiff = curPrecur-prevPrecur;

            int curSNum = Integer.parseInt(str[0]);

            int scanDiff = curSNum-prevSNum;
            long curPos = Long.parseLong(str[1]);
            double mDiff = curPrecur - prevPrecur;

            if( scanDiff==1 || (scanDiff==2 && curDiff==massDiff) )
            //if( massDiff==isoWin || (scanDiff==2 && curDiff==massDiff) )
            {
                posMap.put(curSNum, curPos);
                precursorMap.put(curSNum, curPrecur);
            }
            else if(scanDiff>numIsoWindows || mDiff<0)
            {
                for(double i=prevPrecur+isoWin;i<=lastPrecur;i+=isoWin)
                {
                    prevSNum++;

                    posMap.put(prevSNum, -1);
                    precursorMap.put(prevSNum, i); //prevPrecur+isoWin*(i+1));
                }

                prevSNum++; //skip scan num for the ms1
                for(double i=arr[0];i<curPrecur;i+=isoWin)
                {
                    prevSNum++;

                    posMap.put(prevSNum, -1);
                    precursorMap.put(prevSNum, i); //prevPrecur+isoWin*(i+1));
                }

                posMap.put(curSNum, curPos);
                precursorMap.put(curSNum, curPrecur);
            }
            else
            {
                mDiff = curPrecur - prevPrecur - isoWin;


                    for(int i=0;i<mDiff;i+=isoWin)
                    {
                        prevSNum++;

                        double temp = prevPrecur+isoWin+i;
                        //System.out.println("diff==>>" + mDiff + " " + i + " " + temp + " " + lastPrecur + " "  + isoWin + " " + isoWin*(i+1));
                        //+" " +  prevPrecur+isoWin*(i+1) + indexedFile.getName());
                        if(temp>lastPrecur)
                        {
                            prevSNum++;
                            continue;
                        }

                        posMap.put(prevSNum, -1);
                        precursorMap.put(prevSNum, temp); //prevPrecur+isoWin*(i+1));
                    }

                    posMap.put(curSNum, curPos);
                    precursorMap.put(curSNum, curPrecur); //prevPrecur+isoWin*(i+1));

            }

            prevPrecur = curPrecur;
            prevSNum = curSNum;
        }

        System.out.println(posMap);
        return posMap;
    }
*/

/*
    public MSIndex buildIndex() throws IOException
    {
        String eachLine;

	String temp;

	long pos = file.getFilePointer();

        while( (eachLine=file.readLine()) != null )
        {
                if(eachLine.startsWith("S"))
                {
			temp = eachLine.substring(2);
			temp = temp.substring(0, temp.indexOf("\t"));
                        map.put( Integer.parseInt(temp), pos );
                }
		pos = file.getFilePointer();
        }

        return map;
    }
    */

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

    public TDoubleIntHashMap getRetToScanNumerMap() {
        return retToScanNumerMap;
    }

    public void setRetToScanNumerMap(TDoubleIntHashMap retToScanNumerMap) {
        this.retToScanNumerMap = retToScanNumerMap;
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

    public TDoubleArrayList getRetArr() {
        return retArr;
    }

    public void setRetArr(TDoubleArrayList retArr) {
        this.retArr = retArr;
    }

    public TIntLongHashMap getPosMap() {
        return posMap;
    }

    public void setPosMap(TIntLongHashMap posMap) {
        this.posMap = posMap;
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

    public TIntIntHashMap getIdScanToMs1ScanMap() {
        return idScanToMs1ScanMap;
    }

    public void setIdScanToMs1ScanMap(TIntIntHashMap idScanToMs1ScanMap) {
        this.idScanToMs1ScanMap = idScanToMs1ScanMap;
    }
}

