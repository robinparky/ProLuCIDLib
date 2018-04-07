/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.dia.util;

/**
 *
 * @author rpark
 */


import java.io.*;
import java.util.*;

public class IndexUtil {
    private static final char SPACE = ' ';
    private static int ION_START_INDEX = 3;
    private static final char CARRIAGE_RETURN = '\n';
    private static final char WINDOW_CR = '\r';
    private static final char DOT = '.';


    public static Map<String,IndexedFile> createIndexedHt(String filePath, String extension) throws IOException
    {

        if( !filePath.endsWith(File.separator) )
            filePath += File.separator;

	ArrayList<String> list = FileFilterUtil.getFilesBySuffix(filePath, extension);
        //String[] list = f.list(new RelExFileFilter(extension));

        Map<String, IndexedFile> ht = new Hashtable<>();
        IndexedFile iFile = null;

        String indexFileName;
        File indexFile;

	if(null == list || list.size()<=0) {

		return new Hashtable<String, IndexedFile>();
	}

	if("ms1".equals(extension)) {
		for(Iterator<String> itr=list.iterator(); itr.hasNext(); ) {
		    String each = itr.next();
		    indexFileName = filePath + each + ".index";
		    indexFile = new File(indexFileName);

		    //Create index file
		    if(!indexFile.exists() || indexFile.length()<=0)
		    {
			System.out.println("creating index file " + indexFileName);

			MSIndexFileCreator.createIndexFile(filePath + each); //text.append("Index files are required");
		    }

		    //if index file is corrupted, delete the file and re-try to read it.
		    try {
			iFile = new IndexedFile(indexFile, filePath + each);
		    } catch(FileNotFoundException fnfe) {
                        System.out.println("Error: Spectral files are not found.  If you use mzXML, please use option '-x'");

                    } catch (Exception ioe) {
			System.out.println("re-creating index file " + indexFileName);

			MSIndexFileCreator.createIndexFile(filePath + each); //text.append("Index files are required");
		    }

		    ht.put(each, iFile);
		}

	} else if("mzxml".equals(extension)) {
		for(Iterator<String> itr=list.iterator(); itr.hasNext(); ) {
			String each = itr.next();
		    try {
                        System.out.print("Indexing on " + each + "...");
			iFile = new IndexedFile(filePath + each);
                        System.out.println("done");

		    } catch(FileNotFoundException fnfe) {
                        System.out.println("Error: Spectral file is not found");

                    } catch (Exception ioe) {

			System.out.println("re-creating index file "); // + indexFileName);

			//MSIndexFileCreator.createIndexFile(filePath + each); //text.append("Index files are required");
		    }

                    //System.out.println(ht + " " + iFile + " " + each);
		    ht.put(each, iFile);
		}

	} else if("ms2".equals(extension)) {
		for(Iterator<String> itr=list.iterator(); itr.hasNext(); ) {
			String each = itr.next();
		    indexFileName = filePath + each + ".index";
		    indexFile = new File(indexFileName);

		    if(!indexFile.exists() || indexFile.length()<=0)
		    {
			System.out.println("creating index file " + indexFileName);

			MSIndexFileCreator.createIndexFile(filePath + each); //text.append("Index files are required");
		    }

		    try {
			iFile = new IndexedFile(indexFile, filePath + each);
		    } catch(FileNotFoundException fnfe) {
                        System.out.println("Error: Spectral file is not found");

                    } catch (Exception ioe) {
			ioe.printStackTrace();
			System.out.println(ioe);
			System.out.println("re-creating index file " + indexFileName);

			MSIndexFileCreator.createIndexFile(filePath + each); //text.append("Index files are required");
		    }

		    ht.put(each, iFile);
		}

	}


        return ht;
    }

    public static double[][] getSpectrumArr(IndexedFile iFile, int scanNum) throws IOException
    {
        RandomAccessFile file = iFile.getFile();
        gnu.trove.TIntLongHashMap index = iFile.getMsIndex();
        int[] keys = iFile.getKeys();
        int curIndex = Arrays.binarySearch(keys, scanNum);
        if(curIndex<0) //Cannot find index
            curIndex=-(++curIndex); //Math.abs(++keyIndex);

        BufferedReader br =null;
        String eachLine = null;
        try{
            br = new BufferedReader(new FileReader(iFile.getIndexFile().getPath()));
            while((eachLine=br.readLine()) != null){
            //  System.out.println(scanNum + "\t" + eachLine);
                if(eachLine.startsWith(""+scanNum+"\t")){
                    eachLine=br.readLine();
                    break;
                }
            }

        }catch(Exception e){
            e.printStackTrace();
        }
        if(null == eachLine) {

        }


        if(curIndex>=keys.length)
            curIndex--;

        if(keys.length<=curIndex)
            throw new IOException();

        long startPos = index.get(keys[curIndex]);
        long endPos;

        file.seek(startPos);

        if( (curIndex+1)>=keys.length )
            endPos = file.length();
        else {
          if(null != eachLine) {
            String[] words = eachLine.split("\t");
            endPos = Long.parseLong(words[1]);
          } else {
            endPos = iFile.getFile().length();
          }
        }
            //endPos = index.get(keys[curIndex+1]);

        int byteSize = (int)(endPos-startPos);
        byte[] bytes = new byte[byteSize];
        iFile.getFile().readFully(bytes);

        char ch;
        int pos=0;

        ch = (char)bytes[pos];

        //Remove Z, S, I, D lines
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
        {
            while( ch != CARRIAGE_RETURN )
            {
                pos++;
                ch = (char)bytes[pos];
            }

            pos++;
        }

        StringBuilder mass = new StringBuilder(10);
        StringBuilder intensity = new StringBuilder(10);
        intensity.append('0');

        boolean isMass=true;
        boolean isInt=true;

        int arrSize=0;
        for(int j=pos;j<byteSize;j++)
        {
            if( CARRIAGE_RETURN == (char)bytes[j] )
                arrSize++;
        }

        double[] massArr = new double[arrSize];
        double[] intArr = new double[arrSize];
        double[][] resultArr = new double[2][arrSize];
        resultArr[0] = massArr;
        resultArr[1] = intArr;

        int massIndex=0;

        int spaceCount=0;

        for(int i=pos;i<byteSize;i++)
        {
            ch = (char)bytes[i];
            switch(ch)
            {
                case WINDOW_CR:
                    break;

                case SPACE:
                    spaceCount++;
                    isMass=false;
                    isInt=true;
                    break;

                case CARRIAGE_RETURN:
                    spaceCount=0;
                    isMass=true;
                    isInt=true;

                    //System.out.println(mass.toString() + "\t" + intensity.toString());
                    intArr[massIndex] = Double.parseDouble(intensity.toString());
                    massArr[massIndex++] = Double.parseDouble(mass.toString());

                    mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
                    intensity.delete(0, intensity.length()).append('0');

                    break;

                //case DOT:
                //    isInt=false;

                default:
                    if(spaceCount>=2)
                        break;

                    if(isMass)
                        mass.append(ch);
                    else if(isInt) //remove decimal value of intensity
                        intensity.append(ch);

                    break;
            }
        }


        return resultArr;
    }

    public static void main(String[] args) throws Exception {
      Map<String, IndexedFile> ht = IndexUtil.createIndexedHt("/data/2/rpark/ip2_data/mesingh/NAD_GoGrant_Phospho_HAP/NAD_S12_42_HAPE1_2010_09_17_11_2009/search/projects2017_05_23_12_114353/phospho/temp", "ms2");

      IndexedFile iFile = ht.get("NAD_S12_42_HAPE1_03.ms2");
      System.out.println(ht);
      getSpectrumArr(iFile, 32280);
    }
}
