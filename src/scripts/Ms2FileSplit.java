
package scripts;


import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.Range;

import java.io.*;
import java.util.ArrayList;
import java.util.List;


public class Ms2FileSplit {

    private  int maxMass = 6_000_000;
    private  int minMass = 600_000;
    private StringBuilder cache = new StringBuilder();
    private int firstScan;
    private int lastScan;

    public static void main(String [] args){

        String folder = args[0];

        int numofsplit = Integer.parseInt(args[1]);
        ArrayList<String> ms2files = getMs2Files(folder);
//        new File(new File(folder)+File.separator+"range").mkdir();
       Ms2FileSplit m = new Ms2FileSplit();

	File rangeFolder = new File(folder+File.separator+"range");

	if(rangeFolder.exists()) rangeFolder.delete();
	rangeFolder.mkdir();


        String prefix="";
        if(args.length>2)
            prefix = args[2];
        if(args.length>4)
        {
            int min = Integer.parseInt(args[3]);
            int max = Integer.parseInt(args[4]);
            m.setMinMass(min*1000);
            m.setMaxMass(max*1000);
        }



        for(String ms2File : ms2files){
            m.splitfile(ms2File,numofsplit,6000000, prefix);
        }

                   /*
           deleting null files
           */
        File checkfolder = new File(new File(folder)+File.separator+"range");
        if(checkfolder.exists()) {
            if(checkfolder.isDirectory()) {
                String [] allfiles = checkfolder.list();
                for(String f :allfiles){
                    File check = new File(checkfolder+File.separator+f);
                    long length = check.length();
                    if(length < 100){
                        check.delete();
                    }

                }
            }
        }


    }

    public void splitfile(String ms2file,int splitnum,int maxrange, String prefix){


        File ms2File = new File(ms2file);
        String ms2name = FilenameUtils.removeExtension(ms2File.getName());
        int totalprecmass = this.maxMass;
        List<Integer> range = new ArrayList<>();
        range.add(0);
        for(int i=0;i<totalprecmass;){
            i=i+totalprecmass/splitnum;
            range.add(i);
        }
        BufferedReader br = null;
        try{
          /* File foldercheck = new File(ms2File.getParent()+File.separator+"splitms2");
           if(foldercheck.exists() && foldercheck.isDirectory()){
                for(File f : foldercheck.listFiles()){
                    f.delete();
                }
                foldercheck.delete();
           }*/

            br = new BufferedReader(new FileReader(new File(ms2file)));
            String  eachLine = br.readLine();
            List<BufferedWriter> bout = new ArrayList<>();
            for(int i =0;i<splitnum;i++){
                bout.add(new BufferedWriter(new FileWriter(ms2File.getParent()+File.separator+"range"+File.separator+prefix + ms2name+"_range_"+i+".ms2")));
            }
            for(int i =0;i<splitnum;i++){
               int startRange = range.get(i)-10000;
               if(startRange<0) startRange =0;
               bout.get(i).write("H\tRANGE\t"+(startRange)+"\t"+(range.get(i+1)+10000)+"\n");
            }
            while(eachLine != null){
                if(eachLine.startsWith("S\t")){
                    StringBuffer slinebuffer = new StringBuffer();
                    slinebuffer.append(eachLine+"\n");
                    StringBuffer iLinebuffer = new StringBuffer();
                    StringBuffer massbuffer = new StringBuffer();
                    List<String> zline = new ArrayList<>();
                    String [] words = new String[2];
                    List<Double> precmass = new ArrayList<>();
                    //String [] words = eachLine.split("\t");
                    while((eachLine=br.readLine()) != null && !(eachLine.startsWith("S\t"))){

                        if(eachLine.startsWith("Z\t")){
                            // zlineBuffer.append(eachLine+"\n");
                            zline.add(eachLine+"\n");
                            words = eachLine.split("\t");
                            precmass.add(Double.parseDouble(words[2]));
                        }
                        else if(eachLine.startsWith("I\t")){
                            iLinebuffer.append(eachLine+"\n");
                        }
                        else{
                            massbuffer.append(eachLine+"\n");
                        }
                    }
                    int bufferindex =-1;
                    int zlineindex=0;
                    for(int k=0;k<precmass.size();k++){
                        //(int)((Double.parseDouble(words[3])+0.0005)*1000);

                        int massint = (int)((precmass.get(k)+0.0005)*1000);
                        if(massint<minMass || massint>maxMass) continue;
                        for(int i=0;i<range.size()-1;i++){

                            Range<Integer> r = Range.between(range.get(i), range.get(i+1));
                            if(r.contains(massint)){
                                bufferindex=i;
                                zlineindex=k;
                                break;
                            }
                            else{
                                continue;
                            }
                        }
                        if(bufferindex!=-1) bout.get(bufferindex).write(slinebuffer.toString()+iLinebuffer.toString()+zline.get(zlineindex)+massbuffer);
                    }
                  /* StringBuffer sbtotal = buffer[bufferindex];
                   if(sbtotal == null){
                       sbtotal=sb;
                       buffer[bufferindex]=sbtotal;
                   }else{
                        buffer[bufferindex]=sbtotal.append(sb);
                   }*/
                    continue;
                }

                eachLine=br.readLine();

            }
            for(BufferedWriter out:bout){
                out.close();
            }

           /*new File(ms2File.getParent()+File.separator+"splitms2").mkdir();
           for(int i=0;i<buffer.length;i++){
               if(buffer[i] == null){
                   continue;
               }
               else{
                   File f = new File(ms2File.getParent()+File.separator+"splitms2"+File.separator+ms2name+"_"+i+".ms2");
                   BufferedWriter out = new BufferedWriter(new FileWriter(f));
                   out.write(buffer[i].toString());
                   out.close();
               }
           }*/


        }catch (Exception e){
            e.printStackTrace();
        }
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



    public static String getFileName(String filename)
    {
        int first = filename.lastIndexOf(File.separatorChar);
        int last = filename.lastIndexOf(".");
        return filename.substring(first,last);
    }

    public int getMaxMass() {
        return maxMass;
    }

    public void setMaxMass(int maxMass) {
        this.maxMass = maxMass;
    }

    public int getMinMass() {
        return minMass;
    }

    public void setMinMass(int minMass) {
        this.minMass = minMass;
    }

}

