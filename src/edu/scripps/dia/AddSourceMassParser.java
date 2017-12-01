package edu.scripps.dia;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Titus Jung titusj@scripps.edu on 12/1/17.
 */
public class AddSourceMassParser {


    public static void main(String [] args) throws IOException {
        String sqtFile = args[0];
        String target = args[1];
        String szPath = args[2];
        ParseSource(sqtFile,target,szPath);
    }

    public static void createSimpleSZMap(String szpath, Map<Integer,String> szMap) throws IOException {
        BufferedReader br2 = new BufferedReader(new FileReader(szpath));
        String s,mass ;
        mass = "";
        boolean getMassMode = false;
        String sline;
        String [] arr2=null;
        Map<Integer,String> scanMap;
        Map<Integer,Map<Integer,String>> csMap;
        while((s=br2.readLine())!=null)
        {
            String [] arr = s.split("\t");
            String dataString = arr[0];
            if(dataString.charAt(0)=='Z')
            {
                String name  = arr[0];
                String slineArr[] = arr2;
                int scanNo = Integer.parseInt(slineArr[1]);
                String zlineArr[] = arr;
                int cs = Integer.parseInt(zlineArr[1]);
                mass = zlineArr[2];
                szMap.put(scanNo,mass);
            }
            arr2 = arr;
        }
        br2.close();
    }
    public static void ParseSource(String source, String output, String szfile) throws IOException {
        Map<Integer,String> szMap = new HashMap<>();
        createSimpleSZMap(szfile, szMap);
        BufferedReader br = new BufferedReader(new FileReader(source));
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        String line;
        while((line=br.readLine())!=null)
        {
            String [] arr = line.split("\t");
            int scanNo = Integer.parseInt(arr[1]);
            int cs = Integer.parseInt(arr[4]);
            String [] arr2 = arr[6].split(" ");
            String path = arr2[arr2.length-1];
            String mass =szMap.get(scanNo);
            String [] writeArr = line.split("\t");

            bw.write(line);
            bw.write("\t");
            bw.write(mass);
            bw.newLine();
        }
        bw.close();
        br.close();
    }
}
