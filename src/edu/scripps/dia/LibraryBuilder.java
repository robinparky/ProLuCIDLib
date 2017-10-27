package edu.scripps.dia;

import edu.scripps.dia.models.IndexedMSFile;
import edu.scripps.dia.models.ProcessedPeptide;
import edu.scripps.pms.mspid.ProcessedPeakList;
import edu.scripps.pms.mspid.TerminalModification;
import gnu.trove.TIntIterator;
import org.jdom.JDOMException;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;

/**
 * Created by yateslab on 10/25/17.
 */
public class LibraryBuilder {

    protected List<ProcessedPeptide> peptideList = new ArrayList<>();

    protected Map<String,IndexedMSFile> ms2FileMap = new HashMap<>();


    protected   enum StorageType {SQLITE,PRINT}

    protected StorageType storageType;

    protected String searchXML;
    protected SearchParams sparams;
    protected String terminalString;

    public static void main(String [] args) throws IOException, JDOMException, SQLException {

        String dtaSelectFile = args[0];
        String searchXML = args[1];
        String outputFile = args[2];
        LibraryBuilder builder = new LibraryBuilder(dtaSelectFile,searchXML, StorageType.PRINT);
        builder.printLibraryFile(outputFile);
    }

    public LibraryBuilder()
    {

    }

    public LibraryBuilder(String dtaSelectFile, String searchXML, StorageType storageType) throws IOException, JDOMException {

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

    public void printLibraryFile(String output ) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        int i=0;
        for(ProcessedPeptide peptide: peptideList)
        {

            getSpectrum(peptide,bw);
            if(i%1000==0)
            {
                float percent = i/(float)peptideList.size() *100;
                System.out.println(percent+"% done");
            }
            i++;
        }

        bw.close();
    }

    protected void getSpectrum(ProcessedPeptide peptide, Appendable ap) throws IOException {
        for (String filePath : peptide.getFileScanMap().keySet()) {
            IndexedMSFile file = getMS2File(filePath);
            for (TIntIterator itr = peptide.getFileScanMap().get(filePath).iterator(); itr.hasNext(); ) {
                int scan = itr.next();
                file.goToScan(scan);
                String line;
                boolean start = false;
                readLoop:
                while ((line = file.readLine()) != null) {
                    char firstChar = line.charAt(0);
                    switch (firstChar) {
                        case 'S':
                            if (!start) {
                                start = true;
                                ap.append(line);
                            }
                            else break readLoop;
                            break;
                        default:
                            ap.append(line);
                    }
                    ap.append("\n");
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

}
