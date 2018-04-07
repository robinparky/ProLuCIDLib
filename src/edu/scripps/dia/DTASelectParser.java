package edu.scripps.dia;

import edu.scripps.dia.models.ProcessedPeptide;
import gnu.trove.TIntLongHashMap;

import java.io.*;
import java.util.*;

/**
 * Created by yateslab on 7/18/17.
 */
public class DTASelectParser {


    private Map<String, ProcessedPeptide> processedPeptideMap = new HashMap<>();



    private String dtaselectfilterfilepath;
    private String workDirectory;

    private int seqIndex = -1;
    private int fileNameIndex = -1;

    private int descriptiveIndex = -1;
    private int locusIndex = -1;
    private int theoreticalMassIndex=-1;

    public static void main(String [] args) throws IOException
    {
        DTASelectParser parser = new DTASelectParser(args[0]);
        parser.readDTASelectPepFile();
        Map<String, ProcessedPeptide> map = parser.getProcessedPeptideMap();
        System.out.println("ILDSVGIEADDDR2 has "+map.get("K.ILDSVGIEADDDR.L2").getFileScanSet().size()+ " duplicates");
    }

    public DTASelectParser(String filepath)
    {
        dtaselectfilterfilepath = filepath;
        workDirectory =dtaselectfilterfilepath.substring(0, dtaselectfilterfilepath.lastIndexOf(File.separatorChar));
    }

    public void printProteins(Appendable  out) throws IOException
    {
        for(Iterator<ProcessedPeptide> itr = processedPeptideMap.values().iterator(); itr.hasNext(); )
        {
            ProcessedPeptide peptide = itr.next();
            out.append(peptide.getSequence()+"\n");
            for(Iterator<String> jtr = peptide.proteinIterator(); jtr.hasNext();)
            {
                String proteinInfo = jtr.next();
                out.append("\t"+proteinInfo+"\n");
            }
        }
    }

    public void readDTASelectPepFile(String dtaSelectPepFilePath) throws IOException
    {
        dtaselectfilterfilepath = dtaSelectPepFilePath;
        readDTASelectPepFile();
    }

    public void readDTASelectPepFile() throws IOException
    {
        boolean start = false;
        int peptideRowLength=0;
        int proteinRowLength =0;
        BufferedReader br = new BufferedReader(new FileReader(dtaselectfilterfilepath));
        String line;
        String locus = "";
        String description = "";
        while((line=br.readLine())!=null)
        {
            if(!start && line.startsWith("Locus"))
            {
                String arr [] = line.split("\t");
                proteinRowLength = arr.length;
                processProteinHeader(arr);
            }
            else if(!start && line.startsWith("Unique") )
            {
                String arr [] = line.split("\t");
                peptideRowLength =arr.length;
                start = true;
                processPeptideHeader(arr);
                continue;
            }
            if(start)
            {
                String arr [] = line.split("\t");
                if(arr.length ==peptideRowLength)
                {
                    processPeptideLine(arr,locus,description);
                }
                else if(arr.length==proteinRowLength)
                {
                    locus = arr[locusIndex];
                    description = arr[descriptiveIndex];

                }

            }
        }
        br.close();
    }

    private void processPeptideLine(String [] peptideRow, String locus, String description) throws IOException
    {
        String fileString = peptideRow[fileNameIndex];

        String regex = "\\.";
        String [] fileArr = fileString.split(regex);
        String filePath = workDirectory + File.separatorChar+ fileArr[0];
        int scanNum = Integer.parseInt(fileArr[1]);

        int cs = Integer.parseInt(fileArr[3]);
        float xcorr = Float.parseFloat(peptideRow[2]);
        float deltaCN = Float.parseFloat(peptideRow[3]);
        //String fileScan = filePath+fileArr[1];

        float theorMass= Float.parseFloat(peptideRow[theoreticalMassIndex]);

        String seq =  peptideRow[seqIndex];
        String peptideKey =  seq+ fileArr[3];
        ProcessedPeptide peptide=  processedPeptideMap.get(peptideKey);
        if( peptide !=null)
        {
            peptide.addFileScan(filePath,scanNum);
        }
        else
        {
            peptide = new ProcessedPeptide(seq,theorMass,cs,filePath,scanNum);
            processedPeptideMap.put(peptideKey,peptide);
        }
        peptide.addProtienInfo(locus,description);

    }
    private void processProteinHeader(String [] proteinHeader)
    {
        int i=0;
        for(String str:proteinHeader)
        {
            if(str.contentEquals( "Descriptive Name"))
            {
                descriptiveIndex =i;
            }
            else if(str.contentEquals( "Locus"))
            {
                locusIndex = i;
            }
            i++;
        }
    }
    private void processPeptideHeader(String [] peptideHeader )
    {
        int i=0;
        for(String str:peptideHeader)
        {
            if(str.contentEquals( "Sequence"))
            {
                seqIndex =i;
            }
            else if(str.contentEquals( "FileName"))
            {
                fileNameIndex = i;
            }
            else if(str.contentEquals("CalcM+H+"))
            {
                theoreticalMassIndex = i;
            }
            i++;
        }
    }

    public List<ProcessedPeptide> getProcessedPetideList()
    {
        return new ArrayList<>(processedPeptideMap.values());
    }

    public Map<String, ProcessedPeptide> getProcessedPeptideMap() {
        return processedPeptideMap;
    }

    public String getDtaselectfilterfilepath() {
        return dtaselectfilterfilepath;
    }
}
