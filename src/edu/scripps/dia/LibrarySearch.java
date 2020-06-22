package edu.scripps.dia;

import edu.scripps.dia.util.SpectrumReader;
import edu.scripps.pms.mspid.*;

import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;
import org.jdom.JDOMException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


/**
 * Created by yateslab on 10/12/17.
 */
public class LibrarySearch {
    public static final float DECOY_DIFF = 8;



    public static void rangedIndexSearchMultiThreaded(String ms2path, String paramPath, int numThreads) throws Exception {
        SearchParams sp = new SearchParams(paramPath);
        String libPath = sp.getDatabaseName();
        LibraryIndexer libraryIndexer = new LibraryIndexer(libPath);

        File ms2File = new File(ms2path);
        String ms2DirectoryStr = ms2File.getAbsoluteFile().getParent();

        String msPath = ms2DirectoryStr+File.separator+ms2File.getName();
        String msName = ms2path.substring(0,ms2path.lastIndexOf('.'));
        LibrarySearchEngine clse = new LibrarySearchEngine(ms2path, sp);
        clse.calcRange();
        clse.setLibraryIndexer(libraryIndexer);
        clse.queryRangeFromIndex();
        List<Thread> threadList = new ArrayList<>();
        clse.setOutputPath(msName+".sqt");
        for(int i=0 ; i< numThreads; i++) {
            Thread th = new Thread(new LibrarySearchThread(clse,msName));
            th.start();
            threadList.add(th);
        }
        for(int i=0 ; i< numThreads; i++) {
            threadList.get(i).join();
        }
        clse.clear();
        clse.close();


    }


    public static void main(String[] args) throws Exception {
        //rangedSearched(args);
        //rangedIndexSearch(args[0],args[1]);
        if(args.length>3)
        {
            String paramPath = args[args.length-2];
            int numThreads = Integer.parseInt(args[args.length-1]);
            for(int i=0; i<args.length-2; i++)
            {
                String ms2Path = args[i];
                rangedIndexSearchMultiThreaded(ms2Path,paramPath,numThreads);
            }
        }
        else
        {
            String ms2Path = args[0];
            String paramPath = args[1];
            int numThreads = Integer.parseInt(args[2]);
            rangedIndexSearchMultiThreaded(ms2Path,paramPath,numThreads);
        }

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
          //  sr.setPrelimScoreHits(sphList.subList(0,size));
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
