package edu.scripps.dia.util;

import java.io.*;

/**
 * Created by Titus Jung titusj@scripps.edu on 10/27/17.
 */
public class MSIndexer {


    public static void indexByZline(String ms2file) throws IOException {
        RandomAccessFile raf = new RandomAccessFile(ms2file,"r");
        BufferedWriter bw = new BufferedWriter(new FileWriter(ms2file+".zindex"));
        String line;
        String loLine ="'", z = "" , mass = "", actType = "", retTime ="";
        long pos =0;
        boolean print = false;
        while((line=raf.readLine())!=null)
        {
            char first = line.charAt(0);
            String [] arr;
            switch (first){
                case 'S':
                    arr = line.split("\t");
                    pos = raf.getFilePointer();
                    loLine = arr[1];
                    break;
                case 'Z':
                    arr = line.split("\t");
                     z = arr[1];
                     mass = arr[2];
                     print = true;
                    break;
                case 'I':
                    arr = line.split("\t");
                    String iString = arr[1];
                    if(iString.startsWith("ActivationType"))
                    {
                        actType = arr[2];
                    }
                    else if(iString.startsWith("RetTime"))
                    {
                        retTime = arr[2];
                    }
                    break;
            }
            if(Character.isDigit(first) && print)
            {
               bw.append(loLine).append("\t").append(Long.toString(pos)).append("\t").append(mass).append("\t").
                       append(retTime).append("\t").append(actType).append("\t").append(loLine).append("\t")
                       .append(z);
               bw.newLine();
               print = false;
            }
        }
        bw.close();
    }

    public static void main(String [] args) throws IOException {
        indexByZline(args[0]);
    }


}
