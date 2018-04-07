import edu.scripps.dia.models.IndexedMSFile;
import edu.scripps.dia.models.ProcessedPeptide;
import gnu.trove.TIntIterator;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by Titus Jung titusj@scripps.edu on 1/2/18.
 */
public class ExtractMs2 {

    protected static Map<String,IndexedMSFile> ms2FileMap = new HashMap<>();

    public static void main(String [] args) throws IOException {
           String ms2file = "/home/yateslab/project_data/LibraryBuilder/0301extractms2/QE-C_Ecoli0_5ug_Piercemix_1-2";
           String scanListFile = "/home/yateslab/project_data/LibraryBuilder/0301extractms2/unidentified_pep_from_prolucid_scans.txt";
           String output = "/home/yateslab/project_data/LibraryBuilder/0301extractms2/out.ms2";
           List<Integer> scanList = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(scanListFile));
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        String line;
        while((line=br.readLine())!=null)
        {
            scanList.add(Integer.parseInt(line));
        }
           getSpectrum(scanList,ms2file,bw);
        bw.close();
        br.close();
    }


    public static void getSpectrum(List<Integer> scanList, String filepath, Appendable ap) throws IOException {
        IndexedMSFile file = getMS2File(filepath);
        for (int scan: scanList ) {
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
    private static IndexedMSFile getMS2File(String filepath) throws IOException {
        IndexedMSFile result = ms2FileMap.get(filepath);
        if(result == null)
        {
            result = new IndexedMSFile(filepath);
            ms2FileMap.put(filepath,result);
        }
        return result;
    }
}
