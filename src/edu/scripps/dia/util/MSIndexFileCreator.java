
/*
 * MSIndexFileCreator.java
 *
 * Created on March 25, 2005, 10:36 AM
 */

package edu.scripps.dia.util;

import java.io.*;

public class MSIndexFileCreator
{
  private static final String MS1_FILE = "ms1";
  private static final String MS2_FILE = "ms2";
  private static final String MSZM = "mszm";
  private String filePath;

  private static void printError()
  {
    System.out.println("Usage : java MSIndexFileCreator [-option]");
    System.out.println("-a\t Generate index files for all ms files in the same folder");
    System.out.println("-d\t Generate index files for all ms files in the given folder");
    System.out.println("-f\t ms file name");
  }

  public static void main(String[] args)
    throws IOException
  {
    if ((args.length < 1) || (("-f".equals(args[0])) && (args.length < 2)) || (("-d".equals(args[0])) && (args.length<2)))
    {
      printError();

      return;
    }

    System.out.println("indexing files...");

    if ("-a".equals(args[0]))
    {
      File f = new File(".");
      String[] list = f.list();

      for (int i = 0; i < list.length; ++i)
      {
        if ((null == list[i]) || ((!(list[i].endsWith("ms1"))) && (!(list[i].endsWith("ms2"))) && (!(list[i].endsWith("mszm")))))
          continue;
        createIndexFile(list[i]);
        System.out.println(list[i] + " was indexed.");
      }

      return;
    }
   if("-d".equals(args[0])){
       File f = new File(args[1]);
       String[] list = f.list();

       for (int i = 0; i < list.length; ++i)
       {
           if ((null == list[i]) || ((!(list[i].endsWith("ms1"))) && (!(list[i].endsWith("ms2"))) && (!(list[i].endsWith("mszm")))))
           continue;
           createIndexFile(args[1] + File.separator + list[i]);
           System.out.println(list[i] + " was indexed.");
       }

       return;
   }
    if ("-f".equals(args[0]))
    {
      createIndexFile(args[1]);
      System.out.println(args[1] + " was indexed.");
    }
    else
    {
      printError();
      return;
    }

    System.out.println("indexing completed");
  }

  public MSIndexFileCreator(String filePath)
  {
    this.filePath = filePath;
  }

  public static void createIndexFile(String fileName) throws IOException
  {
    if ((fileName.endsWith("ms1")) || (fileName.endsWith("mszm")))
      createMS1IndexFile(fileName);
    else if (fileName.endsWith("ms2"))
      createMS2IndexFile(fileName);
  }

  private static void createMS1IndexFile(String fileName)
    throws IOException
  {
    File file = new File(fileName);
    InputStream fisIn = new FileInputStream(fileName);

    BufferedReader br = new BufferedReader(new FileReader(fileName));

    FileOutputStream out = new FileOutputStream(fileName + ".index");
    PrintStream p = new PrintStream(out);

    int size = (int)file.length();

    byte[] byteBuffer = new byte[size];

    fisIn.read(byteBuffer, 0, size);
    StringBuffer sb = new StringBuffer();

    for (int i = 0; i < byteBuffer.length; ++i)
    {
      if (((char)byteBuffer[i] != 'S') || ((char)byteBuffer[(i - 1)] != '\n'))
        continue;
      long pos = i;
      int j = i;
      char ch = (char)byteBuffer[j];
      int tabCount = 0;

      while (ch != '\n')
      {
        if (ch == '\t') {
          ++tabCount;
        }
        if (tabCount == 2)
        {
          sb.append(ch);
        }

        ch = (char)byteBuffer[(++j)];
      }

      p.print(Integer.parseInt(sb.toString().trim()));

      p.print("\t");
      p.println(pos);

      sb.delete(0, sb.length());
    }

    fisIn.close();
    p.close();
    out.close();
   br.close();
  }

  private static void createMS2IndexFile(String fileName)
    throws IOException
  {
    File file = new File(fileName);
    InputStream fisIn = new FileInputStream(fileName);

    BufferedReader br = new BufferedReader(new FileReader(fileName));

    FileOutputStream out = new FileOutputStream(fileName + ".index");
    PrintStream p = new PrintStream(out);

    int size = (int)file.length();

    byte[] byteBuffer = new byte[size];

    fisIn.read(byteBuffer, 0, size);
    StringBuffer sb = new StringBuffer();
    StringBuffer precurSb = new StringBuffer();

    for (int i = 0; i < byteBuffer.length; ++i)
    {
      if (((char)byteBuffer[i] != 'S') || ((char)byteBuffer[(i - 1)] != '\n'))
        continue;
      long pos = i;
      int j = i;
      char ch = (char)byteBuffer[j];
      int tabCount = 0;

      while (ch != '\n')
      {
        if (ch == '\t') {
          ++tabCount;
        }
        if (tabCount == 2)
        {
          sb.append(ch);
        }
        else if (tabCount == 3)
        {
          precurSb.append(ch);
        }

        ch = (char)byteBuffer[(++j)];
      }

      p.print(Integer.parseInt(sb.toString().trim()));
      p.print("\t");
      p.print(pos);
      p.print("\t");
      p.println(Float.parseFloat(precurSb.toString().trim()));

      sb.delete(0, sb.length());
      precurSb.delete(0, precurSb.length());
    }

    p.close();
    out.close();
    fisIn.close();
   br.close();
  }

  public static void getSpectrum(String fileName, String scanNum) throws IOException {



//        String spectrumPath = spath + File.separator + "../.." + File.separator + "spectra/" + fileName + ".ms2.index";

   File f = new File(fileName);
   if(!f.exists())
          createIndexFile(fileName);


/*
       BufferedReader br = null;
   try {
       br = new BufferedReader(new FileReader(spectrumPath));

       String eachLine;
       while( (eachLine = br.readLine()) != null) {
           if(eachLine.startsWith(scan))
               break;
       }


       String[] arr = eachLine.split("\t");
       String nextLine = br.readLine();

       RandomAccessFile rfile = new RandomAccessFile(path, "r");
       long startPos = Long.parseLong(arr[1]);
       rfile.seek(startPos);

       int byteSize=-1;

       if(null == nextLine) {
           byteSize = (int)(rfile.length() - startPos);
       } else {
           String[] nextArr = nextLine.split("\t");
           long nextPos = Long.parseLong(nextArr[1]);
           byteSize = (int)(nextPos - startPos);
       }

       byte[] bytes = new byte[byteSize];
       rfile.readFully(bytes);

       //    for(byte b : bytes) {
       //	System.out.print((char)b);
       //  }


       int pos = 0;
       char ch = (char)bytes[pos];
       while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
       {
           while( ch != ProfuseConstants.CARRIAGE_RETURN )
           {
               pos++;
               ch = (char)bytes[pos];
           }

           pos++;
       }

       int arrSize=0;
       for(int j=pos;j<byteSize;j++)
       {
           if( ProfuseConstants.CARRIAGE_RETURN == (char)bytes[j] )
               arrSize++;

       }

       massArr = new double[arrSize];
       intArr = new double[arrSize];


       StringBuilder mass = new StringBuilder(10);
       StringBuilder intensity = new StringBuilder(10);
       intensity.append('0');

       boolean isMass=true;
       boolean isInt=true;

       int massIndex=0;

       for(int j=pos;j<byteSize;j++)
       {
           ch = (char)bytes[j];

           switch(ch)
           {
               case ProfuseConstants.WINDOW_CR:
                   break;

               case ProfuseConstants.SPACE:
                   isMass=false;
                   isInt=true;
                   break;

               case ProfuseConstants.CARRIAGE_RETURN:
                   isMass=true;
                   isInt=false;

                   intArr[massIndex] = Long.parseLong(intensity.toString());
                   massArr[massIndex++] = Double.parseDouble(mass.toString());
                   mass.delete(0, mass.length());
                   intensity.delete(0, intensity.length()).append('0');

                   break;

               case ProfuseConstants.DOT:
                   isInt=false;

               default:
                   if(isMass)
                       mass.append(ch);
                   else if(isInt)
                       intensity.append(ch);

                   break;
           }

       }

   } catch (Exception e) {
       try {   if(null != br) br.close(); }
       catch(IOException ie) { }
   }

       //this.getRequest().setAttribute("massArr", massArr);
       //this.getRequest().setAttribute("intArr", intArr);

*/
  }
}
