package edu.scripps.dia;

import edu.scripps.dia.models.IndexedMSFile;
import edu.scripps.dia.models.ProcessedPeptide;
import edu.scripps.pms.mspid.*;
import gnu.trove.TFloatArrayList;
import gnu.trove.TIntArrayList;
import gnu.trove.TIntIterator;
import org.jdom.JDOMException;

import java.io.*;
import java.util.*;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
/**
 * Created by yateslab on 7/19/17.
 */
public class ConsensusLibraryBuilder  extends LibraryBuilder{



    private int [] spectra = new int[6_000_000];



    public static void main(String [] args) throws IOException, SQLException, JDOMException {

        String dtaSelectFile = args[0];
        String searchXML = args[1];
        String outputFile = args[2];
        ConsensusLibraryBuilder clb = new ConsensusLibraryBuilder(dtaSelectFile, searchXML, StorageType.SQLITE );
        clb.printConsensusFile(outputFile);
      //  clb.uploadConsesusFile(outputFile);
    }

    public ConsensusLibraryBuilder(DTASelectParser parser)
    {
        peptideList = parser.getProcessedPetideList();
    }

    public ConsensusLibraryBuilder(List<ProcessedPeptide> peptideList)
    {
        this.peptideList = peptideList;
    }

    public ConsensusLibraryBuilder(String dtaSelectFile, String searchXML, StorageType storageType ) throws IOException, JDOMException {
        DTASelectParser parser = new DTASelectParser(dtaSelectFile);
        parser.readDTASelectPepFile();
        sparams = new SearchParams(searchXML);
        peptideList = parser.getProcessedPetideList();
        this.storageType = storageType;
        StringBuilder sb = new StringBuilder();
        for(Iterator<TerminalModification> itr = sparams.getCTermDiffMods(); itr.hasNext();)
        {
            TerminalModification mod = itr.next();
            sb.append("+").append(mod.getInfoAsSymbol()).append(":").append(mod.getMassShift());
        }
        for(Iterator<TerminalModification> itr = sparams.getNTermDiffMods(); itr.hasNext();)
        {
            TerminalModification mod = itr.next();
            sb.append("+").append(mod.getInfoAsSymbol()).append(":").append(mod.getMassShift());
        }
        sb.append("+staticC:").append(sparams.getStaticCTermMod());
        sb.append("+staticN:").append(sparams.getStaticNTermMod());

        terminalString = sb.toString();
    }

    /*public void buildConsensuFile(String output) throws IOException
    {
        BufferedWriter bw = new BufferedWriter(new FileWriter(output))
        buildConsensusFile(bw);
        bw.close();
    }*/
    public void uploadConsesusFile(String dbPath) throws SQLException, IOException {
        LibrarySQLiteManager lsm = new LibrarySQLiteManager(dbPath);
        int i=0;
        for(ProcessedPeptide peptide: peptideList)
        {
            buildConsensusSpectrum(peptide);
            TFloatArrayList massList = new TFloatArrayList();
            TIntArrayList intensityList = new TIntArrayList();
            retrievePeaks(massList,intensityList);
            String seq = peptide.getSequence();
            int cs = peptide.getChargeState();
            String seqKey = seq+"."+Integer.toString(cs)+"."+terminalString;
            int massKey= (int)(peptide.getTheoreticalMass()*1000);
            lsm.loadSpectra(seqKey,massKey,cs,massList,intensityList);


            if((i++)%100==0) System.out.print(".");
            //i++;
        }
    }


    public void printConsensusFile(String output ) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        int i=0;
        for(ProcessedPeptide peptide: peptideList)
        {
            bw.append("S");
            bw.newLine();

            bw.append("I\tSeq\t"+peptide.getSequence());
            bw.newLine();

            sparams.getCTermDiffMods();
            sparams.getNTermDiffMods();
            bw.append("I\t").append("TermMods\t").append(terminalString);
            bw.newLine();

            for(String protein: peptide.getProteinInfo())
            {
                bw.append("I\tProtein\t"+protein);
                bw.newLine();
            }

            bw.append("Z\t"+peptide.getChargeState()+"\t"+peptide.getTheoreticalMass());
            bw.newLine();
            buildConsensusSpectrum(peptide);
            printSpectra(bw);
            bw.flush();
            if(i%100==0) System.out.print(".");
            i++;
        }

        bw.close();
    }

    private void printSpectra(Appendable ap) throws IOException
    {
        for(int i=0; i<6_000_000;i++)
        {
            if(spectra[i]!=0)
            {
                float mass = i/1000.0f;
                ap.append(Float.toString(mass)).append(" ").append(Integer.toString(spectra[i])).append("\n");
            }
        }
    }


    private void buildConsensusSpectrum(ProcessedPeptide peptide) throws IOException {

        spectra = new int[6_000_000];
        for(String filePath: peptide.getFileScanMap().keySet())
        {
            IndexedMSFile file = getMS2File(filePath);
            for(TIntIterator itr = peptide.getFileScanMap().get(filePath).iterator(); itr.hasNext();)
            {
                int scan = itr.next();
                file.goToScan(scan);
                String line;
                boolean start = false;
                readLoop: while((line = file.readLine())!=null)
                {
                    char firstChar = line.charAt(0);
                    switch(firstChar)
                    {
                        case 'S':
                            if(!start) start = true;
                            else break readLoop;
                            break;
                        case 'I':

                            break;
                        case 'Z':

                            break;
                        default:
                            if(Character.isDigit(firstChar))
                            {
                                String [] specArray = line.split(" ");
                                double mass = Double.parseDouble(specArray[0]);
                                int intensity =Integer.parseInt(specArray[1]);
                                int massIndex = (int)(mass*1000);
                                spectra[massIndex] += intensity;
                            }
                            break;
                    }
                }
            }
        }
    }

    private IndexedMSFile getMS2File(String filepath) throws IOException {
        IndexedMSFile result = ms2FileMap.get(filepath);
        if(result == null)
        {
            result = new IndexedMSFile(filepath);
            ms2FileMap.put(filepath,result);
        }
        return result;
    }

    public void closeAll()
    {
        for(IndexedMSFile file: ms2FileMap.values())
        {
            try {
                file.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void retrievePeaks(TFloatArrayList massList, TIntArrayList intensityList)
    {
        for(int i=0; i<6_000_000;i++)
        {
            if(spectra[i]!=0)
            {
                float mass = i/1000.0f;
                int intensity = spectra[i];
                massList.add(mass);
                intensityList.add(intensity);
            }
        }
    }


}
