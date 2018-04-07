package edu.scripps.dia.util;

import javax.activation.DataHandler;
import javax.activation.FileDataSource;
import java.io.*;
import java.util.Enumeration;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 *
 * @author Robin Park
 * @version $Id: FileUtil.java,v 1.7 2013/11/26 01:03:57 cvsuser Exp $
 *
 */
public class FileUtil
{
	
	
    public static void main(String[] args) throws Exception {

        //FileUtil.copy( "/home/rpark/testLibCreation.txt", "/home/rpark/test1.txt", false);
        
        //int count = FileUtil.countLines("/data/2/rpark/ip2_data/rpark/Cambridge_Isotope/030714ApoAIWTcmt_2014_04_30_12_24597/spectra/030714ApoAIWTcmt.ms1.index");
        //System.out.println("===" + count);
        
            FileUtil.copy(new File("/home/rampuria/testLibCreation.txt"), new File("/home/rampuria/data/testLibCreation.txt"));
        
        /*FileUtil.replaceCometDBString("/data/2/rpark/ip2_data/rpark/Cambridge_Isotope/030714ApoAIN15CMT_2014_04_30_12_24595/search/comet2014_08_08_17_67605/comet.params", 
                "/data/2/rpark/ip2_data/rpark/Cambridge_Isotope/030714ApoAIN15CMT_2014_04_30_12_24595/search/comet2014_08_08_17_67605/sequest.params",
                "/opt/applications/yates/dbase", "/home/rpark");        
*/
    }

    public FileUtil()
    {
    }

    public static void copy(String in, String out, boolean isAppending) throws IOException {
        copy(new File(in), new File(out), isAppending);
    }

    public static void copy(String in, String out) throws IOException {
        copy(new File(in), new File(out), false);
    }

    public static void copy(File in, File out) throws IOException {
	copy(in, out, false);

    }

    
    public static void copy(File in, File out, boolean isAppending) throws IOException {
        FileInputStream fis  = new FileInputStream(in);
        FileOutputStream fos = new FileOutputStream(out, isAppending);

        byte[] buf = new byte[4096];
        int i = 0;

        while((i=fis.read(buf))!=-1) {
          fos.write(buf, 0, i);
        }

        fis.close();
        fos.close();
    }
    
    public static boolean copyFile(String sourceDir, String fileName,
			String targetDir) throws IOException {
		File source = new File(sourceDir, fileName);
		File target = new File(targetDir, fileName);
		boolean success = FileUtil.streamFileToFile(source, target);
		return success;
	}
    
    public static boolean copyFile(String sourceDir, String fileName,
			String targetDir,String newFileName) throws IOException {
		File source = new File(sourceDir, fileName);
		File target = new File(targetDir, newFileName);
		boolean success = FileUtil.streamFileToFile(source, target);
		return success;
	}
    
    public static boolean fileRename(String sourceDir, String srcName, String newName)
			throws IOException {
		File source = new File(sourceDir, srcName);
		File changedName = new File(sourceDir, newName);
		boolean success = source.renameTo(changedName);
		return success;
	}
    
    public static boolean streamFileToFile(File source, File target)
			throws IOException {
		return handlerStreamFileToFile(source, target);
	}
    
    public static boolean handlerStreamFileToFile(File source, File target)
			throws IOException {
		FileOutputStream stmOut = new FileOutputStream(target);
		try {
			handlerStreamFileToOutputStream(source, stmOut);
		} finally {
			stmOut.close();
		}
		boolean success = target.setLastModified(source.lastModified());
		success = target.setReadOnly();
		return success;
	}
    
    public static void handlerStreamFileToOutputStream(File source,
			OutputStream target) throws IOException {
		FileDataSource dataSource = new FileDataSource(source);
		DataHandler handler = new DataHandler(dataSource);
		handler.writeTo(target);
	}

	public static void expand(File in, File out, boolean isAppending) throws IOException{
		FileInputStream fis = new FileInputStream(in);
		DataInputStream din = new DataInputStream(fis);
		BufferedReader br = new BufferedReader(new InputStreamReader(din));
		int indexIn = in.toString().lastIndexOf("/");
		int indexOu = out.toString().lastIndexOf("/");
		String inPath = in.toString().substring(0,indexIn);
		String ouPath = out.toString().substring(0,indexOu);

		//FileOutputStream fos = new FileOutputStream(out, isAppending);
		FileWriter outFile = new FileWriter(out);
		PrintWriter outW = new PrintWriter(outFile);

		outW.println("#!/bin/bash");
		outW.println("touch " + out.toString() + ".trans\n");
//		fos.write("touch " + out.toString() + ".trans\n");
		String strLine;
		while((strLine = br.readLine()) != null){
			outW.println(strLine);
		}
		outW.println("\nif [ \"$?\" -eq \"0\" ]; then   #chekc it is done");
		outW.println("rm " + out.toString() + ".trans");
                outW.println("touch " +out.toString()+ ".done");
                outW.println("else");
                outW.println("rm " + out.toString() + ".trans");
		outW.println("fi");

/*		
		fos.write("\nif [ \"$?\" -eq \"0\" ]; then   #chekc it is done");
		fos.write("rm " + out.toString() + ".trans");
		fos.write("touch " +out.toString()+ ".done");
		fos.write("else");
		fos.write("rm " + out.toString);
*/

		fis.close();
		din.close();
		outW.close();
	}
	public static int removeAll(String path, String suffix){
		File parent = new File(path);
		int found = 0;
		String[] children = parent.list();
		for(int i=0;i<children.length; i++){
			if(children[i].endsWith(suffix)){
				File toRemove = new File(path+File.separator+children[i]);
				toRemove.delete();
				found++;
			}
		}
		return found;
	}

    public static void createExecutableScript(String scriptname, 
		    String scriptcontent) throws IOException {
        File script = new File(scriptname);
        PrintStream ps = new PrintStream(script);
        ps.print(scriptcontent);
        script.setExecutable(true);
        ps.close();
 
    }

    public static void appendToFile(String filename, String str) throws IOException {
	FileOutputStream fos = new FileOutputStream(filename, true);
	DataOutputStream out   = new DataOutputStream(fos);
	out.writeBytes(str);
	out.flush();
	out.close();
	fos.close();
    }


    //This method is using new lib from jdk and supposed to be faster then old lib.
    //But before using it, make sure it is working correctly.
    //Currently we are not using this method.
    /*
    public static void copy(FileInputStream source, FileOutputStream dest) throws IOException {
         FileChannel in = null, out = null;
         try
         {
              in = source.getChannel();
              out = dest.getChannel();

              long size = in.size();
              MappedByteBuffer buf = in.map(FileChannel.MapMode.READ_ONLY, 0, size);

              out.write(buf);

         }
         catch(IOException e)
         {
             throw new IOException(e.toString());
         }
         finally {
              if (in != null)          in.close();
              if (out != null)     out.close();
         }
    }
    */
    public static boolean deleteDir(String dir) {
	return deleteDir(new File(dir));
    }

    public static boolean deleteDir(File dir) {
	if (dir.isDirectory()) {
	    String[] children = dir.list();
	    for (int i=0; i<children.length; i++) {
		boolean success = deleteDir(new File(dir, children[i]));
		if (!success) {
		    return false;
		}
	    }
	}

	// The directory is now empty so delete it
	return dir.delete();
    }


    public static String lastLine(String fileName)
    {
	//Vector v = tail(fileName, lineCount, 2000);
	Vector v = tail(fileName, 1, 2000);
	return v.get(0).toString();
    }

    public static Vector tail(String fileName, int lineCount)
    {
	//Vector v = tail(fileName, lineCount, 2000);
	return tail(fileName, lineCount, 2000);
    }
    /**
     * Given a byte array this method:
     * a. creates a String out of it
     * b. reverses the string
     * c. extracts the lines
     * d. characters in extracted line will be in reverse order, 
     *    so it reverses the line just before storing in Vector.
     *     
     *  On extracting required numer of lines, this method returns TRUE, 
     *  Else it returns FALSE.
     *   
     * @param bytearray
     * @param lineCount
     * @param lastNlines
     * @return
     */
    private static boolean parseLinesFromLast(byte[] bytearray, int lineCount, Vector lastNlines)
    {
	String lastNChars = new String (bytearray);
	StringBuffer sb = new StringBuffer(lastNChars);
	lastNChars = sb.reverse().toString();
	StringTokenizer tokens= new StringTokenizer(lastNChars,"\n");
	while(tokens.hasMoreTokens())
	{
	    StringBuffer sbLine = new StringBuffer((String)tokens.nextToken());			
	    lastNlines.add(sbLine.reverse().toString());
	    if(lastNlines.size() == lineCount)
	    {
		return true;//indicates we got 'lineCount' lines
	    }
	}
	return false; //indicates didn't read 'lineCount' lines
    }

    /**
     * Reads last N lines from the given file. File reading is done in chunks.
     * 
     * Constraints:
     * 1 Minimize the number of file reads -- Avoid reading the complete file
     * to get last few lines.
     * 2 Minimize the JVM in-memory usage -- Avoid storing the complete file 
     * info in in-memory.
     *
     * Approach: Read a chunk of characters from end of file. One chunk should
     * contain multiple lines. Reverse this chunk and extract the lines. 
     * Repeat this until you get required number of last N lines. In this way 
     * we read and store only the required part of the file.
     * 
     * 1 Create a RandomAccessFile.
     * 2 Get the position of last character using (i.e length-1). Let this be curPos.
     * 3 Move the cursor to fromPos = (curPos - chunkSize). Use seek().
     * 4 If fromPos is less than or equal to ZERO then go to step-5. Else go to step-6
     * 5 Read characters from beginning of file to curPos. Go to step-9.
     * 6 Read 'chunksize' characters from fromPos.
     * 7 Extract the lines. On reading required N lines go to step-9.
     * 8 Repeat step 3 to 7 until 
     *			a. N lines are read.
     *		OR
     *			b. All lines are read when num of lines in file is less than N. 
     * Last line may be a incomplete, so discard it. Modify curPos appropriately.
     * 9 Exit. Got N lines or less than that.
     *
     * @param fileName
     * @param lineCount
     * @param chunkSize
     * @return
     */
    public static Vector tail(String fileName, int lineCount, int chunkSize)
    {
	RandomAccessFile raf = null;
	try
	{
	    raf = new RandomAccessFile(fileName,"r");
	    Vector lastNlines = new Vector();			
	    int delta=0;
	    long curPos = raf.length() - 1;
	    long fromPos;
	    byte[] bytearray;
	    while(true)
	    {				
		fromPos = curPos - chunkSize;
		//System.out.println(curPos);
		//System.out.println(fromPos);				
		if(fromPos <= 0)
		{
		    raf.seek(0);
		    bytearray = new byte[(int)curPos];
		    raf.readFully(bytearray);
		    parseLinesFromLast(bytearray, lineCount, lastNlines);
		    break;
		}
		else
		{					
		    raf.seek(fromPos);
		    bytearray = new byte[chunkSize];
		    raf.readFully(bytearray);
		    if(parseLinesFromLast(bytearray, lineCount, lastNlines))
		    {
			break;
		    }
		    delta = ((String)lastNlines.get(lastNlines.size()-1)).length();
		    lastNlines.remove(lastNlines.size()-1);
		    curPos = fromPos + delta;	
		}
	    }
	    Enumeration e = lastNlines.elements();
//	    while(e.hasMoreElements())
//	    {
//		System.out.println(e.nextElement());
//	    }			
	    return lastNlines;
	}
	catch(Exception e)
	{
	    e.printStackTrace();
	    return null;
	} finally {
		try { if(null != raf) raf.close(); } catch(Exception e) {}
	}
    }	


    public static boolean rename(String target, String dest) {
	File f1 = new File(target);
	File f2 = new File(dest);

	boolean success = f1.renameTo(f2);

	return success;
    }
    
    public static void replaceString(String fname, String oldPattern, String replPattern) {
	String line;
	StringBuffer sb = new StringBuffer();
	try {
            FileInputStream fis = new FileInputStream(fname);
            BufferedReader reader=new BufferedReader ( new InputStreamReader(fis));
            while((line = reader.readLine()) != null) {
                    line = line.replaceAll(oldPattern, replPattern);
                    sb.append(line+"\n");
            }
            reader.close();
            BufferedWriter out=new BufferedWriter ( new FileWriter(fname));
            out.write(sb.toString());
            out.close();
	} catch (Throwable e) {
            System.err.println("Error: " + e);
	}
    }

    public static void replaceSequestDBString(String fname, String oldPattern, String replPattern) {
	String line;
	StringBuffer sb = new StringBuffer();
	try {
            FileInputStream fis = new FileInputStream(fname);
            BufferedReader reader=new BufferedReader ( new InputStreamReader(fis));
            while((line = reader.readLine()) != null) {
            
                if(line.startsWith("database_name"))
                    line = line.replaceAll(oldPattern, replPattern);
                    
                sb.append(line+"\n");
            }
            reader.close();
            BufferedWriter out=new BufferedWriter ( new FileWriter(fname));
            out.write(sb.toString());
            out.close();
	} catch (Throwable e) {
            System.err.println("Error: " + e);
	}
    }
    
    public static void replaceBlazmassDBString(String fname, String outFilename, String oldPattern, String replPattern) {

	String line;
	StringBuffer sb = new StringBuffer();
	try {
            FileInputStream fis = new FileInputStream(fname);
            BufferedReader reader=new BufferedReader ( new InputStreamReader(fis));
            while((line = reader.readLine()) != null) {
            
                if(line.startsWith("[BLAZMASS]")) {
                    line = "[SEQUEST]";
		    sb.append(line+"\n");
		}
			
                if(line.startsWith("database_name")) {
                    line = line.replaceAll(oldPattern, replPattern);
		    sb.append(line+"\n");
		    break;
		}
                    
            }
            reader.close();
            BufferedWriter out=new BufferedWriter ( new FileWriter(outFilename));
            out.write(sb.toString());
            out.close();
	} catch (Throwable e) {
            System.err.println("Error: " + e);
	}
    }

    public static void replaceCometDBString(String fname, String outFilename, String oldPattern, String replPattern) {

	String line;
	StringBuffer sb = new StringBuffer();
        sb.append("[SEQUEST]\n");
        
	try {
            FileInputStream fis = new FileInputStream(fname);
            BufferedReader reader=new BufferedReader ( new InputStreamReader(fis));
            while((line = reader.readLine()) != null) {
            
                if(line.startsWith("database_name")) {
                    line = line.replaceAll(oldPattern, replPattern);
		    sb.append(line+"\n");
		    break;
		}
                    
            }
            reader.close();
            BufferedWriter out=new BufferedWriter ( new FileWriter(outFilename));
            out.write(sb.toString());
            out.close();
	} catch (Throwable e) {
            System.err.println("Error: " + e);
	}
    }
    
    public static int countLines(String filename)  {
        InputStream is = null;
        try {
            is = new BufferedInputStream(new FileInputStream(filename));
            byte[] c = new byte[1024];
            int count = 0;
            int readChars = 0;
            boolean empty = true;
            while ((readChars = is.read(c)) != -1) {
                empty = false;
                for (int i = 0; i < readChars; ++i) {
                    if (c[i] == '\n') {
                        ++count;
                    }
                }
            }
            return (count == 0 && !empty) ? 1 : count;
        } catch (Exception e) {
            e.printStackTrace();
        }finally {
            try {
                if(null == is)
                    is.close();
            } catch(Exception e) {                
            }
        }
        
        return 0;
    }
    
    public static void makeDir(String targetDir) {
    	try{
    		File dataFile = new File(targetDir);
    		if (!dataFile.exists()) {
    			dataFile.mkdirs();
    		}
    	}catch (Exception ex) {
			System.out.println("makeDir :: "+ex.getMessage());
		}

	}
    
    public static void writeJSON(String jsonData, String jsonFileName) {

		try {
			BufferedWriter writer = null;
			File f = new File(jsonFileName);
			if (f.exists()) {
				f.delete();
			}
			writer = new BufferedWriter(new FileWriter(jsonFileName, true));
//			writer.write(JsonWriter.formatJson(jsonData));
			writer.write(jsonData);
			writer.write("\n");
			writer.close();
		} catch (IOException ex) {
			System.out.println("writeJSON :: "+ex.getMessage());
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
    
    public static void appendFile(File srcFile, File dstFile) {

		FileWriter fstream = null;
		BufferedWriter out = null;

		try {

			fstream = new FileWriter(srcFile, true);
			out = new BufferedWriter(fstream);

		} catch (IOException e1) {

			e1.printStackTrace();
		}

		FileInputStream fis;

		try {
			
			fis = new FileInputStream(dstFile);
			
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));

			String aLine;
			
			while ((aLine = in.readLine()) != null) {
				out.write(aLine);
				out.newLine();
			}
			
			in.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
