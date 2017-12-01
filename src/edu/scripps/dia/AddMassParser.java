package edu.scripps.dia;

import gnu.trove.TIntObjectHashMap;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Titus Jung titusj@scripps.edu on 11/30/17.
 */
public class AddMassParser {

    public static void main(String [] args) throws IOException {
        String sqtFile = args[0];
        String target = args[1];
        String szPath = args[2];
        ParseAll(sqtFile,target,szPath);
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
            String [] arr = s.split(":");
            String dataString = arr[1];
            if(dataString.charAt(0)=='Z')
            {
                String name  = arr[0];
                String slineArr[] = arr2[1].split("\t");
                int scanNo = Integer.parseInt(slineArr[1]);
                String zlineArr[] = arr[1].split("\t");
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
            int scanNo = Integer.parseInt(arr[2]);
            int cs = Integer.parseInt(arr[4]);
            String [] arr2 = arr[6].split(" ");
            String path = arr2[arr2.length-1];
            String mass =szMap.get(scanNo);
            bw.write(line);
            bw.write("\t");
            bw.write(mass);
            bw.newLine();
        }
        bw.close();
        br.close();
    }




    public static void createSZMap(String szpath, Map<String,Map<Integer,Map<Integer,String>>> szMap ) throws IOException {
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
            String [] arr = s.split(":");
            String dataString = arr[1];
            if(dataString.charAt(0)=='Z')
            {
                String name  = arr[0];
                String slineArr[] = arr2[1].split("\t");
                int scanNo = Integer.parseInt(slineArr[1]);
                String zlineArr[] = arr[1].split("\t");
                int cs = Integer.parseInt(zlineArr[1]);
                 mass = zlineArr[2];

                csMap = szMap.get(name);
                if(csMap==null)
                {
                    scanMap = new HashMap<>();
                    scanMap.put(scanNo,mass);
                    csMap =new HashMap<>();
                    csMap.put(cs,scanMap);
                    szMap.put(name,csMap);
                }
                else
                {
                    scanMap = csMap.get(cs);
                    if(scanMap == null)
                    {
                        scanMap = new HashMap<>();
                        scanMap.put(scanNo,mass);
                        csMap.put(cs,scanMap);
                    }
                    else
                    {
                        scanMap.put(scanNo,mass);
                    }
                }

            }
            arr2 = arr;
        }
        br2.close();

    }

    public static String getMass(Map<String,Map<Integer,Map<Integer,String>>> szMap,String path, int cs, int scanNo )
    {
        String name = path.substring(path.lastIndexOf('/')+1,path.length());
        return szMap.get(name).get(cs).get(scanNo);
    }

    public static void ParseAll(String txtAllFile, String output, String sz) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(txtAllFile));
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        Map<String,Map<Integer,Map<Integer,String>>> szMap = new HashMap<>();
        createSZMap(sz,szMap);
        String line;
        while((line=br.readLine())!=null)
        {
            String [] arr = line.split("\t");
            int scanNo = Integer.parseInt(arr[2]);
            int cs = Integer.parseInt(arr[4]);
            String [] arr2 = arr[6].split(" ");
            String path = arr2[arr2.length-1];
            String mass =getMass(szMap,path,cs,scanNo);
            bw.write(line);
            bw.write("\t");
            bw.write(mass);
            bw.newLine();
        }
        bw.close();
        br.close();

        /*
        String line;
        while((line=br.readLine())!=null)
        {
            String [] arr = line.split("\t");
            int scanNo = Integer.parseInt(arr[2]);
            int cs = Integer.parseInt(arr[4]);
            String [] arr2 = arr[6].split(" ");
            String path = arr2[arr2.length-1];
            String mass =getMass(scanNo,path,cs,sz);
            bw.write(line);
            bw.write("\t");
            bw.write(mass);
            bw.newLine();
        }
        bw.close();
        br.close();*/
    }




    public static void ParseSQT(String sqtFile, String target) throws IOException {
   /*     BufferedReader br = new BufferedReader(new FileReader(sqtFile));
        BufferedWriter bw = new BufferedWriter(new FileWriter(target));
        String line;
        while((line=br.readLine())!=null)
        {
            if(line.startsWith("M"))
            {
                String [] arr = line.split("\t");
                String scanNo = arr[2];
                int scanNum = Integer.parseInt(scanNo);
                String [] arr2 = arr[arr.length-1].split(" ");
                String path = arr2[arr2.length-1];
                int insert = arr.length-1;
                if(arr[arr.length-2].equals("Decoy"))
                {
                    insert --;
                }
                String mass = getMass(scanNum,path, -1);
                for(int i=0; i<arr.length; i++)
                {
                    if(i==insert)
                    {
                        bw.write(mass);
                        bw.write("\t");
                    }
                    bw.write(arr[i]);
                    bw.write("\t");
                }
                bw.newLine();

            }
            else
            {
                bw.write(line);
                bw.newLine();
            }
        }
        bw.close();
        br.close();*/

    }

    public static String getMass(int scanNo, String path, int cs, String szPath) throws IOException {
        String name = path.substring(path.lastIndexOf('/')+1,path.length());
        BufferedReader br2 = new BufferedReader(new FileReader(szPath));
        String s,mass ;
        mass = "";
        boolean getMassMode = false;
        while((s=br2.readLine())!=null)
        {
            String [] arr0 = s.split(":");
            String dataString = arr0[1];
            if(name.equals(arr0[0]))
            {
                if(getMassMode && dataString.charAt(0)=='Z')
                {
                    String [] arr = s.split("\t");
                    int z = Integer.parseInt(arr[1]);
                    if(z!=cs)
                    {
                        getMassMode = false;
                        continue;
                    }
                    mass = arr[2];
                    break;
                }
                if(dataString.charAt(0)=='S')
                {
                    String [] arr = s.split("\t");
                    String code = arr[1];
                    int icode = Integer.parseInt(code);
                    if(icode==scanNo)
                    {
                        getMassMode = true;
                        continue;
                    }
                }
            }
        }
        br2.close();
        return mass;
    }

}
