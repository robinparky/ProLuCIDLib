package edu.scripps.dia;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import edu.scripps.dia.util.SpectrumReader;
import edu.scripps.pms.mspid.*;
import edu.scripps.pms.mspid.SearchParams;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;
import org.jdom.JDOMException;
import org.omg.PortableServer.SERVANT_RETENTION_POLICY_ID;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;


/**
 * Created by yateslab on 10/12/17.
 */
public class LibrarySearch {
    public static final double DECOY_DIFF = 50.0;

    public static void simpleSearch(String [] args) throws Exception {
        String ms2Path = args[0];
        String paramsPath = args[1];

        LibrarySearchEngine clse = new LibrarySearchEngine(ms2Path, paramsPath);
        clse.calcRange();
        SearchParams params = clse.getSearchParams();
        String libraryPath = params.getDatabaseName();
        String path = ms2Path.substring(0,ms2Path.lastIndexOf("."));
        String output = path+"_"+libraryPath.substring(0,libraryPath.lastIndexOf("."))+".sqt";
        //clse.readMS2TargetFile(libraryPath);
        clse.readMs2TargetDirectory(libraryPath);
      //  List<SearchResult> results = new ArrayList<>();
        List<PeakList> peakLists = clse.getPeakLists();
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));

        for(PeakList peakList: peakLists)
        {
            for (Iterator<Zline> it = peakList.getZlines(); it.hasNext(); ) {
                Zline zline = it.next();
                SearchResult r = clse.search(peakList,zline);
                bw.write(r.outputResults());
                bw.newLine();
            }
        }
        bw.close();
    }

    public static void rangedSearched(String [] args) throws Exception {
        String ms2Path = args[0];
        String paramsPath = args[1];
        SearchParams sp = new SearchParams(paramsPath);
        String libraryPath = sp.getDatabaseName();
       // String output = args[3];
        File libDirectory= new File(libraryPath);
        File ms2Directory = new File(ms2Path);
        String[] ms2Array = ms2Directory.list();
        String ms2DirectoryStr = ms2Path;

        if(ms2Directory.isFile())
        {
            ms2DirectoryStr = ms2Directory.getAbsoluteFile().getParent();
            ms2Array = new String[]{ms2Directory.getName()};
        }
        List<SearchResult> srList = new ArrayList<>();

        String[] libArray = libDirectory.list();
        for(String ms2Files: ms2Array)
        {
            if(ms2Files.endsWith("ms2"))
            {
                String msPath = ms2DirectoryStr+File.separator+ms2Files;
                String msName = ms2Files.substring(0,ms2Files.lastIndexOf('.'));
                LibrarySearchEngine clse = new LibrarySearchEngine(msPath, paramsPath);
                clse.calcRange();
                clse.getLibraryRanges(libArray,libraryPath);
\
                List<String> targetFiles = clse.getFilesInRange();


                String output = ms2DirectoryStr+File.separator+msName+".sqt";

                BufferedWriter bw = new BufferedWriter(new FileWriter(output));
                bw.write(SearchResult.SLINEHEADER);
                bw.newLine();
                bw.write(SearchResult.MLINEHEADER);
                bw.newLine();
                for(String libFiles: targetFiles)
                {


                    if(libFiles.endsWith("ms2"))
                    {
                        //int index = libFiles.lastIndexOf(".");

                       // System.out.println(output);
                        String libPath = libraryPath+File.separator+libFiles;
                        clse.readMS2TargetFile(libPath);
                        System.out.println("searching: "+libFiles);
                        List<PeakList> peakLists = clse.getPeakLists();

                        for(PeakList peakList: peakLists) {
                            for (Iterator<Zline> it = peakList.getZlines(); it.hasNext(); ) {
                                Zline zline = it.next();
                                SearchResult r = clse.search(peakList, zline);
                                r.setMs2Name(msName);
                                srList.add(r);
                                if (srList.size() % 100 == 0) {
                                    for (SearchResult sr : srList) {
                                        bw.write(sr.outputResults());
                                    }
                                    srList = new ArrayList<>();
                                }
                                //bw.newLine();
                            }
                        }
                        for (SearchResult sr : srList) {
                            bw.write(sr.outputResults());
                        }
                        srList = new ArrayList<>();

                    }
                    clse.clear();

                }

                bw.close();
                System.out.println("written to: "+output);
            }

        }
    }


    public static void main(String[] args) throws Exception {
        rangedSearched(args);
         //simpleSearch(args);
        //inverseSearch(args);
    }

    public static void inverseSearch(String[] args) throws IOException, JDOMException {
        DistributionCalculator dc = new DistributionCalculator();
        String srcDirectoryPath = args[0];
        String paramsPath = args[1];
        String targetDirectoryPath = args[2];
        String output = args[3];

        edu.scripps.dia.SearchParams params = new edu.scripps.dia.SearchParams(paramsPath);
        Map<String,ProcessedPeakList> pplMap = new HashMap<>();
        Map<String,PriorityQueue<ScoredPeptideHit>> srMap = new HashMap<>();
        MassCalculator mc = new MassCalculator(params);

        File srcDirectory = new File(srcDirectoryPath);
        String[] srcFiles = srcDirectory.list();
        for(String f:srcFiles)
        {
            if(f.endsWith("ms2"))
            {
                SpectrumReader spectrumReader = new SpectrumReader(srcDirectoryPath+File.separatorChar+f,"ms2");
                List<PeakList> peakLists = spectrumReader.getSpectraList();
                for(PeakList pl: peakLists)
                {
                    pl.setFilename(f);
                    ProcessedPeakList ppl = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc);
                    pplMap.put(f+"_"+pl.getHiscan(),ppl);
                }
                spectrumReader.closeDataFile();
            }
        }

        File targetDirectory = new File(targetDirectoryPath);
        String [] targetFiles = targetDirectory.list();

        for(String f:targetFiles)
        {
            if(f.endsWith("ms2"))
            {
                System.out.println(f);
                SpectrumReader spectrumReader = new SpectrumReader(targetDirectoryPath+File.separatorChar+f,"ms2");
                Iterator<PeakList> plIterator= spectrumReader.getSpectra();
                while(plIterator.hasNext())
                {
                    PeakList pl = plIterator.next();
                    pl.setFilename(f);
                    Zline zline = pl.getZlines().next();
                    ProcessedPeakList targetppl = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc);
                    for(Map.Entry<String,ProcessedPeakList> entry: pplMap.entrySet())
                    {
                        ProcessedPeakList sourceppl = entry.getValue();
                        String key = entry.getKey();
                        double sourceMass = sourceppl.getZline().getM2z();


                        if(MassRangeFinder.withinRange(sourceMass,zline.getM2z(),params))
                        {
                            System.out.println(f+"\t"+zline.getM2z());

                            ScoredPeptideHit sph = search(sourceppl,targetppl);
                            if(srMap.containsKey(key))
                            {
                                PriorityQueue<ScoredPeptideHit> srQueue =srMap.get(key);
                                if(srQueue.peek().getPScore()>sph.getPScore() || srQueue.size()<params.getCandidatePeptideThreshold())
                                {
                                    srQueue.add(sph);
                                    if(srQueue.size()>params.getCandidatePeptideThreshold())
                                    {
                                        srQueue.poll();
                                    }
                                }
                            }
                            else
                            {
                                PriorityQueue<ScoredPeptideHit> srQueue = new PriorityQueue<>(new ScoredPeptideHitComparator());
                                srQueue.add(sph);
                                srMap.put(key,srQueue);
                            }

                        }
                    }
                }
                spectrumReader.closeDataFile();

            }
        }

        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        for(Map.Entry<String,PriorityQueue<ScoredPeptideHit>> entry:srMap.entrySet())
        {
            String filename = entry.getKey();
            PriorityQueue queue = entry.getValue();
            ProcessedPeakList ppl = pplMap.get(filename);
            SearchResult sr = new SearchResult(ppl);
            List<ScoredPeptideHit> sphList = new ArrayList<>(queue);
            int size = sphList.size()<SearchResult.getNUMFINALRESULT()?sphList.size():SearchResult.getNUMFINALRESULT();
            sr.setPrelimScoreHits(sphList.subList(0,size));
            bw.append(sr.outputResults());
          //  bw.newLine();
            bw.flush();
        }
        bw.close();
    }

    public static ScoredPeptideHit search(ProcessedPeakList source, ProcessedPeakList target)  {


        PeakList pl = target.getPeakList();
        ScoredPeptideHit  sph= source.prelimScoreCorrelation(target);
        sph.setMs2CompareValues(pl.getHiscan(),pl.getLoscan(), pl.getFilename());

        return sph;
    }
}
