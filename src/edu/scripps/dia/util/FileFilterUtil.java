package edu.scripps.dia.util;


import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 *
 * @author Robin Park
 * @author Tao Xu
 * @version $Id: FileFilterUtil.java,v 1.10 2013/05/03 21:22:17 cvsuser Exp $
 */
public class FileFilterUtil
{
    public static FileFilter getDirectoryFilter()
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return file.isDirectory();
            }
        };
    }

    public static FileFilter getFileFilter()
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return !file.isDirectory();
                // file.isFile() does not display link files
//                return file.isFile();
            }
        };
    }

    public static FileFilter getSQTFilter()
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return file.getName().endsWith("sqt");
            }
        };
    }

    public static FileFilter getExactFileFilter(final String fileName)
    {
        return new FileFilter()
        {
            public boolean accept(File file)
            {
                return file.getName().equals(fileName);
            }
        };
    }
    // get files by suffix and prefix
    public static ArrayList<String> getFilesBySuffix(String dir, String prefix, String suffix) {
         ArrayList<String> fFiles = new ArrayList<String>();
         File currentDir = new File(dir);
         String [] files = currentDir.list();
	 if(files==null) return fFiles;

       for(int i = 0; i < files.length; i++) {
//                 System.out.println("got file: " + s);
             String s = files[i];
             if(s.endsWith(suffix) && s.startsWith(prefix)) {
                 fFiles.add(s);
             }
         }
         return fFiles;
    }
    // get files by suffix
    public static ArrayList<String> getFilesBySuffix(String dir, String suffix) {
         //System.out.println("getFilesBySuffix in " + dir + " for suffix " + suffix);
         ArrayList<String> suffixFiles = new ArrayList<String>();
         File currentDir = new File(dir);
         String [] files = currentDir.list();
if(files==null) return suffixFiles;
       for(int i = 0; i < files.length; i++) {
//                 System.out.println("got file: " + s);
             String s = files[i];
             if(s.endsWith(suffix)) {
                 suffixFiles.add(s);
             }
         }
         return suffixFiles;
    }

    public static ArrayList<String> getFilesByPrefix(String dir, String prefix) {
         ArrayList<String> prefixFiles = new ArrayList<String>();
         File currentDir = new File(dir);
         String [] files = currentDir.list();
	 if(files==null) return prefixFiles;

       for(int i = 0; i < files.length; i++) {
//                 System.out.println("got file: " + s);
             String s = files[i];
             if(s.startsWith(prefix)) {
                 prefixFiles.add(s);
             }
         }
         return prefixFiles;
    }

    public static ArrayList<String> getFilesBySuffixWithPath(String dir, String suffix) {

         ArrayList<String> suffixFiles = new ArrayList<String>();
         File currentDir = new File(dir);

  //      System.out.println("==" + dir);

         for(String s : currentDir.list()) {
//                 System.out.println("got file: " + s);
             if(s.endsWith(suffix)) {
                 suffixFiles.add(dir + File.separator + s);
             }
         }
         return suffixFiles;
    }

    public static ArrayList<String> getMs2Files(String f) {
        ArrayList<String> ms2files = new ArrayList<String>();
        File myfile = new File(f);
        if(myfile.exists()) {
            if(myfile.isDirectory()) {
                String [] allfiles = myfile.list();
                for(int i = 0; i < allfiles.length; i++) {
                    String file = allfiles[i];
                    if(file != null) {
                        if(file.endsWith(".ms2") || file.endsWith(".mzXML")) {
                            String ms2filename = f + File.separator + file;
                            ms2files.add(ms2filename);
                            //System.out.println("adding " + ms2filename);
                        }
                    }
                }
            } else {
                ms2files.add(f);
            }
        }
        return ms2files;
    }

    // files are sorted from small to big in size
    private static ArrayList<String> getMs2FilesSortedBySize(String f) {
        ArrayList<String> ms2files = new ArrayList<String>();
        ArrayList<SortableFile> sortedms2files = new ArrayList();

        File myfile = new File(f);
        if(myfile.exists()) {
            if(myfile.isDirectory()) {
                String [] allfiles = myfile.list();
                for(int i = 0; i < allfiles.length; i++) {
                    String file = allfiles[i];
                    if(file != null) {
                        if(file.endsWith(".ms2") || file.endsWith(".mzXML")) {
                            String ms2filename = f + File.separator + file;
                            sortedms2files.add(new SortableFile(ms2filename));
                            //System.out.println("adding " + ms2filename);
                        }
                    }
                }
            } else {
                sortedms2files.add(new SortableFile(f));
            }
        }
        Collections.sort(sortedms2files);

        for(Iterator<SortableFile> it = sortedms2files.iterator(); it.hasNext();) {
            SortableFile sf = it.next();

            //System.out.println(sf.getFileName() + " "  + sf.getFileSize() + "\n");
            ms2files.add(sf.getFileName());
        }
        return ms2files;
    }

    // return two list of ms2 files, the first list of small files should be submitted together, while
    // the big ms2 files in the second list should be submitted individually
    // each group of ms2 files is guaranteed to have no more than maxnumspectra or only one file
    public static ArrayList<ArrayList<String>> sortMs2FilesBySize(String dir, int maxnumspectra) {

         ArrayList<ArrayList<String>> groupedms2files = new ArrayList();
         ArrayList<String> ms2files = new ArrayList();
         int numspectraadded = 0;
         //for(Iterator<String> it = getMs2Files(dir).iterator(); it.hasNext();) {
         for(Iterator<String> it = getMs2FilesSortedBySize(dir).iterator(); it.hasNext();) {
             String f = it.next();
             int numspectra = (int) new File(f).length()/5000; // 5000 is the average size of a ms2 spectra

             int newnumspectra = numspectraadded + numspectra;

             if(newnumspectra <= maxnumspectra) {
                 ms2files.add(f);
                 numspectraadded = newnumspectra;
             } else {
                 if(ms2files.size() > 0) {
                     groupedms2files.add(ms2files);
                 }
                 ms2files = new ArrayList();
                 ms2files.add(f);
                 numspectraadded = numspectra;
             }
         }

         // deal with the last group of ms2 files
         if(ms2files.size() > 0) {
             groupedms2files.add(ms2files);
         }
         return groupedms2files;
    }

    public static int countFilesBySuffix(String dir, String suffix) {
         int numfiles = 0;
         File currentDir = new File(dir);
      if(!currentDir.exists()) return numfiles;

         for(String s : currentDir.list()) {
                 //System.out.println("got file: " + s);
             if(s.endsWith(suffix)) {
                 numfiles++;
             }
         }
         return numfiles;
    }

    public static int countFilesBySuffix(String dir, String suffix, String excludeString) {
        int numfiles = 0;
        File currentDir = new File(dir);

        for(String s : currentDir.list()) {
            //System.out.println("got file: " + s);
            if(s.endsWith(suffix) && !s.contains(excludeString)) {
                numfiles++;
            }
        }
        return numfiles;
    }


    public static String getFilePrefix(String filename, String suffix) {
        return filename.indexOf(suffix) == -1? filename : filename.split(suffix)[0];
    }
    public static void copyfiles(String sourceFolder, String targetFolder,
                 String suffix) throws IOException, InterruptedException {

	File f = new File(targetFolder);
	if(!f.exists())
	    f.mkdir();

        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
        Iterator<String> fileit = files.iterator();
        //System.out.println("===copy " + files.size() + " *" + suffix + " files from " + sourceFolder + " to " + targetFolder);
        while(fileit.hasNext()) {
            String file = fileit.next();
            String sourcefile = sourceFolder + "/" + file;
            String targetfile = targetFolder + "/" + file;
            String lncommand = "cp " + sourcefile + " " + targetfile;
            //System.out.println("===command==" + lncommand);
            Process p = Runtime.getRuntime().exec(lncommand);
            p.waitFor();
	    //System.out.println("===waitfor return value " + p.waitFor());

        }
    }

    public static boolean isHeavyFile(String file, String path) {

	    if(!file.startsWith("H"))
		return false;


	    String lightfile = path + "/" + file.substring(1, file.length());

	    File lf = new File(lightfile);

	    if(lf.exists())
		    return true;
	    else
		    return false;

    }

    public static boolean isMediumFile(String file, String path) {

	    if(!file.startsWith("M"))
		return false;


	    String lightfile = path + "/" + file.substring(1, file.length());

	    File lf = new File(lightfile);

	    if(lf.exists())
		    return true;
	    else
		    return false;

    }


    public static void copyfilesWithoutHeavySqt(String sourceFolder, String targetFolder,
                 String suffix) throws IOException, InterruptedException {

	File f = new File(targetFolder);
	if(!f.exists())
	    f.mkdir();

        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
        Iterator<String> fileit = files.iterator();

	//boolean heavySqtExcluded = false;
        //System.out.println("===copy " + files.size() + " *" + suffix + " files from " + sourceFolder + " to " + targetFolder);
        while(fileit.hasNext()) {
            String file = fileit.next();
            String sourcefile = sourceFolder + "/" + file;

	    //check if it is heavy search
	    if(file.startsWith("H")) {
		String lightfile = sourceFolder + "/" + file.substring(1, file.length());

		File lf = new File(lightfile);

		if(lf.exists()) {
		   // heavySqtExcluded = true;
		    continue;
		}
	    }


            String targetfile = targetFolder + "/" + file;
            String lncommand = "cp " + sourcefile + " " + targetfile;
            //System.out.println("===command==" + lncommand);
            Process p = Runtime.getRuntime().exec(lncommand);
            p.waitFor();
	    //System.out.println("===waitfor return value " + p.waitFor());

        }

	//return heavySqtExcluded;
    }

    public static void copyfilesHeavySqtOnly(String sourceFolder, String targetFolder,
                 String suffix) throws IOException, InterruptedException {

	File f = new File(targetFolder);
	if(!f.exists())
	    f.mkdir();

        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
        Iterator<String> fileit = files.iterator();

	//boolean heavySqtExcluded = false;
        //System.out.println("===copy " + files.size() + " *" + suffix + " files from " + sourceFolder + " to " + targetFolder);
        while(fileit.hasNext()) {
            String file = fileit.next();
            String sourcefile = sourceFolder + "/" + file;

	    //check if it is heavy search
	    if(file.startsWith("H")) {
		String lightfile = sourceFolder + "/" + file.substring(1, file.length());

	//	System.out.println("heavy " + sourcefile);
	//	System.out.println("lheavy " + lightfile);

		File lf = new File(lightfile);

		if(!lf.exists()) {
		   // heavySqtExcluded = true;
		    continue;
		}


                String targetfile = targetFolder + "/" + file;
                String lncommand = "cp " + sourcefile + " " + targetfile;
                //System.out.println("===command==" + lncommand);
                Process p = Runtime.getRuntime().exec(lncommand);
                p.waitFor();
	    }


	    //System.out.println("===waitfor return value " + p.waitFor());

        }

	//return heavySqtExcluded;
    }


    public static void makeHardcopy(String sourceFolder, String targetFolder,
                 String suffix) throws IOException, InterruptedException {

	File f = new File(targetFolder);
	if(!f.exists())
	    f.mkdir();

        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
        Iterator<String> fileit = files.iterator();
        while(fileit.hasNext()) {
            String file = fileit.next();
            String sourcefile = sourceFolder + "/" + file;
            String targetfile = targetFolder + "/" + file;
            //String lncommand = "ln " + "-s " + sourcefile + " " + targetfile;
            //Process p = Runtime.getRuntime().exec(lncommand);
	    //p.waitFor();
	    FileUtil.copy(sourcefile, targetfile);
        }
    }

    public static void makeRelativeSymbolicLinks(String sourceFolder, String targetFolder,
                 String suffix, String relativePath) throws IOException, InterruptedException {
/*
	File f = new File(targetFolder);
	if(!f.exists())
	    f.mkdir();

        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
        Iterator<String> fileit = files.iterator();
        while(fileit.hasNext()) {
            String file = fileit.next();
            String sourcefile = sourceFolder + "/" + file;
            String targetfile = targetFolder + "/" + file;
            String lncommand = "ln " + "-s " +  relativePath + File.separator + file + " " + targetfile;
            //String lncommand = "ln " + "-s " + sourcefile + " " + targetfile;
        //    System.out.println("==command==" + lncommand);
            Process p = Runtime.getRuntime().exec(lncommand);
	    p.waitFor();

        }
*/
	makeRelativeSymbolicLinks(sourceFolder, targetFolder, suffix, relativePath, "");
    }

    public static void makeRelativeSymbolicLinks(String sourceFolder, String targetFolder,
                 String suffix, String relativePath, String prefixtoadd) throws IOException, InterruptedException {
//System.out.println("sourceFolder " + sourceFolder + " targetFolder " + targetFolder + " suffix " + suffix + " relativePath " + relativePath);
	File f = new File(targetFolder);
	if(!f.exists())
	    f.mkdir();

        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
//System.out.println("sourceFolder: " + sourceFolder + " number of " + files.size() + " file s with suffix: " + suffix + " with relativePath " + relativePath);
        Iterator<String> fileit = files.iterator();
        while(fileit.hasNext()) {
            String file = fileit.next();
            //String sourcefile = sourceFolder + "/" + file;
            String targetfile = targetFolder + "/" + prefixtoadd + file;
            String lncommand = "ln " + "-s " +  relativePath + File.separator + file + " " + targetfile;
            //String lncommand = "ln " + "-s " + sourcefile + " " + targetfile;
            //System.out.println("==command==" + lncommand);
            Process p = Runtime.getRuntime().exec(lncommand);
	    p.waitFor();
//System.out.println("Exit value of " + lncommand + " command is " + p.exitValue());
        }
    }

    public static void makeSymbolicLink(String sourceFolder, String targetFolder, String sourceFile, String targetFile) throws IOException, InterruptedException {

	    File f = new File(targetFolder);
	    if(!f.exists())
		    f.mkdir();

	    String sourcefile = sourceFolder + "/" + sourceFile;
	    String targetfile = targetFolder + "/" + targetFile;
	    String lncommand = "ln " + "-s " + sourcefile + " " + targetfile;
	    //    System.out.println("==command==" + lncommand);
	    Process p = Runtime.getRuntime().exec(lncommand);
	    p.waitFor();


    }
    public static void makeSymbolicLink(String sourceFolder, String targetFolder, String file) throws IOException, InterruptedException {

	    File f = new File(targetFolder);
	    if(!f.exists())
		    f.mkdir();

	    String sourcefile = sourceFolder + "/" + file;
	    String targetfile = targetFolder + "/" + file;
	    String lncommand = "ln " + "-s " + sourcefile + " " + targetfile;
	    //    System.out.println("==command==" + lncommand);
	    Process p = Runtime.getRuntime().exec(lncommand);
	    p.waitFor();


    }

    public static void makeSymbolicLinks(String sourceFolder, String targetFolder,
                 String suffix) throws IOException, InterruptedException {

	File f = new File(targetFolder);
	if(!f.exists())
	    f.mkdir();

        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
        Iterator<String> fileit = files.iterator();
        while(fileit.hasNext()) {
            String file = fileit.next();
            String sourcefile = sourceFolder + "/" + file;
            String targetfile = targetFolder + "/" + file;
	    String lncommand = "ln " + "-s " + sourcefile + " " + targetfile;
        //    System.out.println("==command==" + lncommand);
            Process p = Runtime.getRuntime().exec(lncommand);
	    p.waitFor();

        }
    }
    // for heavy search, add a prefix to the name of the links,
    // e.g., link 1.ms2 to prefixtoadd1.ms2
    public static void makeSymbolicLinks(String sourceFolder, String targetFolder,
                 String suffix, String prefixtoadd) throws IOException, InterruptedException {
        ArrayList<String> files = getFilesBySuffix(sourceFolder, suffix);
        Iterator<String> fileit = files.iterator();
        while(fileit.hasNext()) {
            String file = fileit.next();
            String sourcefile = sourceFolder + "/" + file;
            String targetfile = targetFolder + "/" + prefixtoadd + file;
            String lncommand = "ln " + "-s " + sourcefile + " " + targetfile;
            //System.out.println(lncommand);
            Process p = Runtime.getRuntime().exec(lncommand);
	    p.waitFor();
        }
    }

    public static boolean exists(String file) {
        return new File(file).exists();
    }
    private static long estimateNumSpectra(ArrayList<String> files) {
        long numspectra = 0;
        for(Iterator<String> it = files.iterator(); it.hasNext();) {
            numspectra += new File(it.next()).length()/5000;
        }
        return numspectra;
    }
    public static void main(String[] args) throws Exception {


        int num = FileFilterUtil.countFilesBySuffix("/home/rpark/testLibCreation", ".java", "_my_");
        System.out.println(num);

        if(true) return;

        System.out.println("== " + FileFilterUtil.isHeavyFile("H0730-L2-10.ms2", "/data/2/rpark/ip2_data/rpark/Eddie_Liao/Control_L2_2014_09_29_15_27521/search/projects2014_09_29_22_70553") );
    //makeRelativeSymbolicLinks("/home/rpark/test2", "/home/rpark/test2/a", "ms2", "..", "H");

	    ArrayList<ArrayList<String>> alist = FileFilterUtil.sortMs2FilesBySize(args[0], Integer.parseInt(args[1]));
	    for(Iterator<ArrayList<String>> itr = alist.iterator(); itr.hasNext(); ) {
                    ArrayList<String> files = itr.next();
		    System.out.println(files + " " + estimateNumSpectra(files));
		    System.out.println("====");
	    }

    }
}
