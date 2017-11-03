package edu.scripps.dia;

import edu.scripps.dia.util.SpectrumReader;
import edu.scripps.pms.mspid.*;
import edu.scripps.pms.mspid.SearchParams;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;
import org.jdom.JDOMException;

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
        String libraryPath = args[2];
        String output = args[3];
        LibrarySearchEngine clse = new LibrarySearchEngine(ms2Path, paramsPath);
        clse.calcRange();
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
        String libraryPath = args[2];
       // String output = args[3];
        File libDirectory= new File(libraryPath);
        File ms2Directory = new File(ms2Path);


        String[] libArray = libDirectory.list();
        String[] ms2Array = ms2Directory.list();
        for(String ms2Files: ms2Array)
        {
            if(ms2Files.endsWith("ms2"))
            {
                String msPath = ms2Directory+File.separator+ms2Files;
                String msName = ms2Files.substring(0,ms2Files.lastIndexOf('.'));
                LibrarySearchEngine clse = new LibrarySearchEngine(msPath, paramsPath);
                clse.calcRange();

                for(String libFiles: libArray)
                {
                    int index = libFiles.lastIndexOf(".");
                    String output = ms2Path+File.separator+msName+"_"+libFiles.substring(0,index)+".sqt";
                    String libPath = libraryPath+File.separator+libFiles;

                    if(libFiles.endsWith("ms2"))
                    {
                        clse.readMS2TargetFile(libPath);
                        System.out.println("searching: "+libFiles);
                        List<PeakList> peakLists = clse.getPeakLists();
                        BufferedWriter bw = new BufferedWriter(new FileWriter(output));

                        for(PeakList peakList: peakLists)
                        {
                            for (Iterator<Zline> it = peakList.getZlines(); it.hasNext(); ) {
                                Zline zline = it.next();
                                SearchResult r = clse.search(peakList,zline);
                                r.setMs2Name(msName);
                                bw.write(r.outputResults());
                                //bw.newLine();
                            }
                        }
                        bw.close();
                        System.out.println("written to: "+output);
                        clse.clear();
                    }
                }
            }

        }
    }


    public static void main(String[] args) throws Exception {
        rangedSearched(args);
        // simpleSearch(args);
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
                    ProcessedPeakList ppl = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc,true);
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
                    ProcessedPeakList targetppl = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc,true);
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
