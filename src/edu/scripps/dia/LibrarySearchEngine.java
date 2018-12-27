package edu.scripps.dia;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import edu.scripps.dia.util.MzxmlPeakList;
import edu.scripps.dia.util.MzxmlSpectrumReader;
import edu.scripps.dia.util.SpectrumReader;
import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;
import gnu.trove.TIntArrayList;
import gnu.trove.TIntHashSet;
import org.jdom.JDOMException;

import java.io.*;
import java.sql.SQLException;
import java.util.*;

/**
 * Created by yateslab on 8/8/17.
 */
public class LibrarySearchEngine {
    private static int NUMFINALRESULT = 5;
    private RangeMap<Integer,List<String>> massRangeFileMap =  TreeRangeMap.create();
    private int tempMaxCS = 20;
    private List<PeakList> peakLists;
    private List<PeakList> libraryPeakLists = new ArrayList<>();
    private List<ProcessedPeakList>[] libraryPeakListTable = new List[tempMaxCS];
    private List<LibrarySpectra>[] indexedSpectraListTable = new List[tempMaxCS];

    private int [][] massIndex = new int[tempMaxCS][];
    private List<String> [] fileIndex = new List[6_500_000];
    private LibraryIndexer libraryIndexer;
    private SearchParams params;
    private String ms2FilePath;
    //private LibrarySQLiteManager library;
    private MassCalculator mc;
    private int startRange=-1;
    private int endRange=-1;
    private boolean useSQLite = false;
    private static DistributionCalculator dc;
    private int searchCount =0;
    private Set<ProcessedPeakList> cache = new HashSet<>();
    private PeakCache pplCache = new PeakCache();
    private Iterator<PeakList> peakListIterator=null;
    private BufferedWriter bw = null;
    private String outputPath;

    public static void main(String[] args) throws Exception {
        String ms2Path = args[0];
        String paramsPath = args[1];
        String libraryPath = args[2];
        String output = args[3];
        LibrarySearchEngine clse = new LibrarySearchEngine(ms2Path, paramsPath);
        clse.calcRange();
        //clse.readMS2TargetFile(libraryPath);
        clse.readMs2TargetDirectory(libraryPath);
        List<SearchResult> results = new ArrayList<>();
        List<PeakList> peakLists = clse.getPeakLists();
        for(PeakList peakList: peakLists)
        {
            for (Iterator<Zline> it = peakList.getZlines(); it.hasNext(); ) {
                Zline zline = it.next();
                results.add(clse.search(peakList,zline));
            }
        }

        BufferedWriter bw = new BufferedWriter(new FileWriter(output));

        for(SearchResult r: results)
        {
            bw.write(r.outputResults());
        }
        bw.close();
    }
    public LibrarySearchEngine( String searchParamsPath) throws IOException, JDOMException, SQLException {
        params = new SearchParams(searchParamsPath);
        mc = new MassCalculator(params);
        if(dc == null) {
            dc = new DistributionCalculator();
        }
    }

    public LibrarySearchEngine(String ms2Path, String searchParamsPath) throws IOException, JDOMException, SQLException {
        params = new SearchParams(searchParamsPath);
        ms2FilePath = ms2Path;
        mc = new MassCalculator(params);
        if(dc == null) {
            dc = new DistributionCalculator();
        }
    }


    public LibrarySearchEngine(String ms2Path, SearchParams params) {
        this.params = params;
        ms2FilePath = ms2Path;
    }

    private void init() throws Exception {
        peakLists = getSpectra(ms2FilePath);
    }
    public SearchResult searchIndexedDB(PeakList peakList, Zline zline) throws SQLException {
        int numPeaks = peakList.numPeaks();
        searchCount ++;
        PriorityQueue<ScoredPeptideHit> sphQueue = new PriorityQueue<>(new ScoredPeptideHitComparator());
        TIntHashSet set = new TIntHashSet();
        double worseScore = Double.MAX_VALUE;
        int z= zline.getChargeState();
        ProcessedPeakList ppl = new ProcessedPeakList(peakList, zline, params, mc);

        SearchResult searchResult = new SearchResult(ppl);
        // System.out.println(peakList.getLoscan());
        if (numPeaks > params.getMinNumSpectra() && numPeaks < params.getMaxNumSpectra()) {
            float prcMass = (float) zline.getM2z();
            ScoredPeptideHit sph;
            if (prcMass > params.getMinPrecursorMass() && prcMass < params.getMaxPrecursorMass()) {
                TIntArrayList highLimits = new TIntArrayList();
                TIntArrayList lowLimits = new TIntArrayList();
                MassRangeFinder.findRange(prcMass, params, lowLimits, highLimits);
                ppl.preprocess(params.getPreprocess());
                searchResult = new SearchResult(ppl);
                if(massIndex[z]==null) return  searchResult;
                set.clear();
                for (int i = 0; i < highLimits.size(); i++) {
                    int high = highLimits.get(i);
                    int low = lowLimits.get(i);
                    int massLow = massIndex[zline.getChargeState()][low-startRange];
                    int massHigh = massIndex[zline.getChargeState()][high-startRange];
                    //  System.out.println("low "+massLow+" high "+massHigh);

                    for(int j=massLow; j<massHigh; j++)
                    {
                        if(set.contains(j)){
                            continue;
                        }
                        set.add(j);
                        LibrarySpectra spectra = indexedSpectraListTable[zline.getChargeState()].get(j);
                     //   ProcessedPeakList ppl2 =indexedSpectraListTable[zline.getChargeState()].get(j);
                        //ProcessedPeakList ppl2 = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc,true);
                       // PeakList pl = ppl2.getPeakList();
                        sph= ppl.prelimScoreCorrelation(spectra);
                        //cache.add(ppl2);
                        //pplCache.put(ppl.getID(),ppl);

                     //   sph.setMs2CompareValues(pl.getHiscan(),pl.getLoscan(), pl.getFilename());
                        if(sph.getPScore()<worseScore || sphQueue.size()<params.getCandidatePeptideThreshold())
                        {
                            sphQueue.add(sph);
                            if(sphQueue.size()>params.getCandidatePeptideThreshold())
                            {
                                sphQueue.poll();
                            }
                            worseScore = sphQueue.peek().getPrimaryScore();
                        }
                    }

                }

                highLimits = new TIntArrayList();
                lowLimits = new TIntArrayList();
                MassRangeFinder.findRange(prcMass+LibrarySearch.DECOY_DIFF, params, lowLimits, highLimits);
                for (int i = 0; i < highLimits.size(); i++) {
                    int high = highLimits.get(i);
                    int low = lowLimits.get(i);
                    if(low<startRange || high>endRange ) continue;
                    int massLow = massIndex[zline.getChargeState()][low-startRange];
                    int massHigh = massIndex[zline.getChargeState()][high-startRange];
                    // System.out.println("low "+massLow+" high "+massHigh);

                    for(int j=massLow; j<massHigh; j++)
                    {
                        if(set.contains(j)){
                            continue;
                        }
                        set.add(j);
                        set.add(j);
                        LibrarySpectra spectra = indexedSpectraListTable[zline.getChargeState()].get(j);
                        //   ProcessedPeakList ppl2 =indexedSpectraListTable[zline.getChargeState()].get(j);
                        //ProcessedPeakList ppl2 = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc,true);
                        // PeakList pl = ppl2.getPeakList();
                        sph= ppl.prelimScoreCorrelation(spectra);
                        //cache.add(ppl2);
                       // pplCache.put(ppl.getID(),ppl);

                        //   sph.setMs2CompareValues(pl.getHiscan(),pl.getLoscan(), pl.getFilename());
                        if(sph.getPScore()<worseScore || sphQueue.size()<params.getCandidatePeptideThreshold())
                        {
                            sphQueue.add(sph);
                            if(sphQueue.size()>params.getCandidatePeptideThreshold())
                            {
                                sphQueue.poll();
                            }
                            worseScore = sphQueue.peek().getPrimaryScore();
                        }
                        sph.setIsDecoy(true);

                    }

                }

                //ppl.dumpBoolMass();


                List<ScoredPeptideHit> sphList = new ArrayList<>(sphQueue);
                int size = sphList.size()<NUMFINALRESULT?sphList.size():NUMFINALRESULT;
                searchResult.setPrelimScoreHits(sphList,size);
                searchResult.calcScores();
            }
        }
        // System.out.println();
        //searchResult.setFinalResultsForPrelimScores();
        //System.out.println(">>>>cache size is "+cache.size());

        return searchResult;
    }


    public SearchResult search(PeakList peakList, Zline zline) throws SQLException {
        int numPeaks = peakList.numPeaks();
        searchCount ++;
        PriorityQueue<ScoredPeptideHit> sphQueue = new PriorityQueue<>(new ScoredPeptideHitComparator());
        TIntHashSet set = new TIntHashSet();
        double worseScore = Double.MAX_VALUE;
        int z= zline.getChargeState();
        ProcessedPeakList ppl = new ProcessedPeakList(peakList, zline, params, mc);

        SearchResult searchResult = new SearchResult(ppl);
       // System.out.println(peakList.getLoscan());
        if (numPeaks > params.getMinNumSpectra() && numPeaks < params.getMaxNumSpectra()) {
            float prcMass = (float) zline.getM2z();
            ScoredPeptideHit sph;
            if (prcMass > params.getMinPrecursorMass() && prcMass < params.getMaxPrecursorMass()) {
                TIntArrayList highLimits = new TIntArrayList();
                TIntArrayList lowLimits = new TIntArrayList();
                MassRangeFinder.findRange(prcMass, params, lowLimits, highLimits);
                ppl.preprocess(params.getPreprocess());
                searchResult = new SearchResult(ppl);
                if(massIndex[z]==null) return  searchResult;
                set.clear();
                for (int i = 0; i < highLimits.size(); i++) {
                    int high = highLimits.get(i);
                    int low = lowLimits.get(i);
                    int massLow = massIndex[zline.getChargeState()][low-startRange];
                    int massHigh = massIndex[zline.getChargeState()][high-startRange];
                  //  System.out.println("low "+massLow+" high "+massHigh);

                    for(int j=massLow; j<massHigh; j++)
                    {
                        if(set.contains(j)){
                            continue;
                        }
                        set.add(j);
                        ProcessedPeakList ppl2 =libraryPeakListTable[zline.getChargeState()].get(j);
                        //ProcessedPeakList ppl2 = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc,true);
                        PeakList pl = ppl2.getPeakList();
                        sph= ppl.prelimScoreCorrelation(ppl2);
                        //cache.add(ppl2);
                        pplCache.put(ppl.getID(),ppl);
                      //  System.out.println("<<<size " +ppl2.getMassSet().size());
                        // sph = ppl.prelimScoreCorrelation(pl);
                      //  System.out.println(pl.getLoscan()+"\t"+sph.getPScore());
                         sph.setMs2CompareValues(pl.getHiscan(),pl.getLoscan(), pl.getFilename());
                         if(sph.getPScore()<worseScore || sphQueue.size()<params.getCandidatePeptideThreshold())
                         {
                             sphQueue.add(sph);
                             if(sphQueue.size()>params.getCandidatePeptideThreshold())
                             {
                                   sphQueue.poll();
                             }
                             worseScore = sphQueue.peek().getPrimaryScore();
                         }
                    }

                }
                 highLimits = new TIntArrayList();
                 lowLimits = new TIntArrayList();
                MassRangeFinder.findRange(prcMass+LibrarySearch.DECOY_DIFF, params, lowLimits, highLimits);
                for (int i = 0; i < highLimits.size(); i++) {
                    int high = highLimits.get(i);
                    int low = lowLimits.get(i);
                    if(low<startRange || high>endRange ) continue;
                    int massLow = massIndex[zline.getChargeState()][low-startRange];
                    int massHigh = massIndex[zline.getChargeState()][high-startRange];
                    // System.out.println("low "+massLow+" high "+massHigh);

                    for(int j=massLow; j<massHigh; j++)
                    {
                        if(set.contains(j)){
                            continue;
                        }
                        set.add(j);
                        ProcessedPeakList ppl2 =libraryPeakListTable[zline.getChargeState()].get(j);
                        //ProcessedPeakList ppl2 = new ProcessedPeakList(pl,pl.getZlines().next(),params,mc,true);
                        PeakList pl = ppl2.getPeakList();
                        sph= ppl.prelimScoreCorrelation(ppl2);
                        sph.setIsDecoy(true);
                        // sph = ppl.prelimScoreCorrelation(pl);
                        //  System.out.println(pl.getLoscan()+"\t"+sph.getPScore());
                        sph.setMs2CompareValues(pl.getHiscan(),pl.getLoscan(), pl.getFilename());
                        if(sph.getPScore()<worseScore || sphQueue.size()<params.getCandidatePeptideThreshold())
                        {
                            sphQueue.add(sph);
                            if(sphQueue.size()>params.getCandidatePeptideThreshold())
                            {
                                sphQueue.poll();
                            }
                            worseScore = sphQueue.peek().getPrimaryScore();
                        }
                        //cache.add(ppl2);
                        pplCache.put(ppl.getID(),ppl);

                    }

                }

                ppl.dumpBoolMass();


                List<ScoredPeptideHit> sphList = new ArrayList<>(sphQueue);
                int size = sphList.size()<NUMFINALRESULT?sphList.size():NUMFINALRESULT;
                searchResult.setPrelimScoreHits(sphList,size);
            }
        }
       // System.out.println();
        //searchResult.setFinalResultsForPrelimScores();
        //System.out.println(">>>>cache size is "+cache.size());

        return searchResult;
    }

    public ScoredPeptideHit search(ProcessedPeakList source, ProcessedPeakList target) throws SQLException {


        PeakList pl = target.getPeakList();
        ScoredPeptideHit  sph= source.prelimScoreCorrelation(target);
        sph.setMs2CompareValues(pl.getHiscan(),pl.getLoscan(), pl.getFilename());

        return sph;
    }


    public SearchResult inverseSearch(PeakList peakList, Zline zline, String ms2FilePath) throws IOException, SQLException {
        clear();
        readMS2TargetFile(ms2FilePath);
        return search(peakList,zline);
    }


    private float simpleScore(float[] consensusSpectra, float[] massCorr, float[] massSum, int[] peaksInfo) {
        float score = 0;
        int numPeaks = 0;
        int matchedPeaks = 0;
        for (int i = 0; i < consensusSpectra.length; i++) {
            if (consensusSpectra[i] > 0) {
                numPeaks++;
                if (massCorr[i] != 0) {
                    matchedPeaks++;
                    score += consensusSpectra[i] * massCorr[i];
                }
            }
        }
        peaksInfo[0] = numPeaks;
        peaksInfo[1] = matchedPeaks;
        return score;
    }


    public static List<PeakList> getSpectra(String file) throws IOException, JDOMException, Exception {
        ArrayList<PeakList> peaklists = null;
        if (file.endsWith(".mzXML")) {
            MzxmlSpectrumReader sr = new MzxmlSpectrumReader(file);
            ArrayList<PeakList> spectra = new ArrayList<PeakList>(20000);
            for (Iterator<MzxmlPeakList> it = sr.getSpectra(2); it.hasNext(); ) {
                spectra.add(it.next());
            }
            peaklists = spectra;
            sr.closeDataFile();
            sr = null;
        } else {
            SpectrumReader sr = new SpectrumReader(file, "ms2");
            peaklists = sr.getSpectraList();
            sr.closeDataFile();
            sr = null;
        }

        return peaklists;
    }

    public List<PeakList> getPeakLists() throws Exception {
        if(peakLists==null)
            peakLists = getSpectra(ms2FilePath);
        return peakLists;
    }



    public synchronized PeakList getPeakList() throws Exception {
        if(peakLists==null)
            peakLists = getSpectra(ms2FilePath);
        if(peakListIterator==null)
        {
            peakListIterator = peakLists.listIterator();
        }
        if(!peakListIterator.hasNext()) return null;
        return peakListIterator.next();
    }

    public void readMs2TargetDirectory(String ms2DirectoryPath) throws IOException {
        File directory= new File(ms2DirectoryPath);
        String[] fileArray = directory.list();
        for(String fileName: fileArray)
        {
            if(fileName.endsWith("ms2"))
            {
                readMS2TargetFile(ms2DirectoryPath+File.separatorChar+ fileName);
            }
        }

    }

    public void readMS2TargetFile(String ms2FilePath) throws IOException {
        SpectrumReader reader = new SpectrumReader(ms2FilePath,"ms2");
       // List<PeakList> plList = reader.getSpectraList();
        for(Iterator<PeakList> pItr = reader.getSpectra(); pItr.hasNext(); )
        {
            PeakList peakList = pItr.next();

           int z=  peakList.getFirstChargeState();
           peakList.setFilename(ms2FilePath);
           double mass = peakList.getZlines().next().getM2z();
           int massloc = (int)(mass*1000);
           int massLocation = massloc -startRange;
           if(massloc>=endRange || massloc < startRange) continue;
           if(libraryPeakListTable[z]==null)
           {
            //   System.out.println(">>> cs "+z);
               libraryPeakListTable[z] = new ArrayList<>();
               massIndex[z] = new int[endRange - startRange];
           }
           ProcessedPeakList ppl = new ProcessedPeakList(peakList,peakList.getZlines().next(),params,mc);
           libraryPeakListTable[z].add(ppl);
       //    System.out.println(">>> "+z+"\t"+massLocation+"\t"+mass);
           massIndex[z][massLocation]++;
        }
        reader.closeDataFile();

        fillIndex();
        sortTargetSpectra();

    }

    public void readMS2LibraryFile(String ms2FilePath) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(ms2FilePath));
        String line;
        boolean readSpectraMode = false;
        PeakList pl = new PeakList();
        double mass;
        int cs;
        String seq,def, def2;
        while((line=br.readLine())!=null)
        {
            if(line.charAt(0) == 'Z')
            {
                String [] arr = line.split("\t");
                 mass = Double.parseDouble(arr[2]);
                 cs = Integer.parseInt(arr[1]);
                pl = new PeakList();
                readSpectraMode = true;
            }
            else if(!readSpectraMode && line.charAt(0)=='I')
            {
                String [] arr = line.split("\t");
                if(arr[1].equals("Seq"))
                {
                   seq = arr[2];
                }
                else if(arr[1].equals("Protein"))
                {
                     def = arr[2];
                     def2 = arr[3];
                }
                //else if()
            }
            else if(line.charAt(0) == 'S')
            {
                libraryPeakLists.add(pl);
                readSpectraMode = false;
            }
            else if(readSpectraMode &&  Character.isDigit( line.charAt(0)))
            {
                Peak peak = new Peak(line);
                pl.addPeak(peak);
            }
        }
        fillIndex();
        sortTargetSpectra();
    }

    public int getIndex(double mass)
    {
        return (int)(mass*1000)-startRange;
    }

    public void clear()
    {
        for(int j=0; j<tempMaxCS; j++)
        {
            massIndex[j]=null;
            libraryPeakListTable[j] =null;
        }

    }

    public void fillIndex()
    {
        for(int j=0; j<tempMaxCS; j++)
        {
            if(massIndex[j]!=null)
            {
                for(int i=1; i <massIndex[j].length; i++)
                {
                    massIndex[j][i]+= massIndex[j][i-1];
                }
            }

        }

    }

    public void sortLibrarySpectra()
    {
        Comparator<LibrarySpectra> spectraComparator = new Comparator<LibrarySpectra>() {
            @Override
            public int compare(LibrarySpectra spectra, LibrarySpectra spectra2) {
                return Double.compare(spectra.mz,spectra2.mz);
            }
        };
        for(int i=0; i<indexedSpectraListTable.length; i++)
        {
            if(indexedSpectraListTable[i]!=null)
            {
                Collections.sort(indexedSpectraListTable[i],spectraComparator);
            }
        }
    }


    public void sortTargetSpectra()
    {
        Comparator<ProcessedPeakList> peakListComparator = new Comparator<ProcessedPeakList>() {
            @Override
            public int compare(ProcessedPeakList peakList, ProcessedPeakList t1) {
                return Double.compare(peakList.getPrecursorMass(), t1.getPrecursorMass());
            }
        };
        for(int i=0; i<libraryPeakListTable.length; i++)
        {
            if(libraryPeakListTable[i]!=null)
            {
                Collections.sort(libraryPeakListTable[i],peakListComparator);
            }
        }
    }

    public void calcRange() throws IOException {
        String eachLine;
        int startrange=600_000;
        int endrange=6_000_000;
        BufferedReader br = new BufferedReader(new FileReader(ms2FilePath));
        while ((eachLine = br.readLine()).startsWith("H\tRANGE")) {
            String[] words = eachLine.split("\t");
            startrange = Integer.parseInt(words[2]);
            endrange = Integer.parseInt(words[3]);

            if (startrange < 590000) {
                continue;
            }


            if (params.getStaticNTermMod() > 0) {
                startrange = (int) (startrange - params.getStaticNTermMod() * 1000);
            }
            if (params.getStaticCTermMod() > 0) {
                startrange = (int) (startrange - params.getStaticCTermMod() * 1000);
            }
            Iterator itr2 = params.getNTermDiffMods();
            while (itr2.hasNext()) {
                TerminalModification t = (TerminalModification) itr2.next();
                if (t.getMassShift() > 0) {
                    startrange = (int) (startrange - t.getMassShift() * 1000);
                }
            }
            Iterator itr3 = params.getCTermDiffMods();
            while (itr3.hasNext()) {
                TerminalModification t = (TerminalModification) itr3.next();
                if (t.getMassShift() > 0) {
                    startrange = (int) (startrange - t.getMassShift() * 1000);
                }
            }

          //  System.out.println("" + ms2FilePath);
            List<Double> massshift = new ArrayList<>();
            Iterator itr = params.getDiffMods();
            while (itr.hasNext()) {
                DiffMod d = (DiffMod) itr.next();
                massshift.add(d.getMassShift());
            }

            int massshiftno = params.getMaxAlter();
            if (!(massshift.size() == 0)) {
                double maxshift = Collections.max(massshift);
                double minshift = Collections.min(massshift);
                startrange = startrange - (int) ((maxshift) * 1000 * massshiftno);
                if (minshift < 0) {
                    endrange = endrange - (int) (minshift * 1000 * massshiftno);
                }
            }
        }
        this.startRange = startrange;
        this.endRange = endrange+1000;
        br.close();
    }

    public boolean queryRangeFromIndex() throws IOException, SQLException {
        if(endRange==-1 || startRange==-1) calcRange();
        if(libraryIndexer==null) return false;
        System.out.println("Start "+startRange+" end "+endRange);

        List<LibrarySpectra> spectraList = libraryIndexer.querySpectra(startRange,endRange);
        for(LibrarySpectra spectra : spectraList)
        {
            float mass = spectra.mz;

            int massloc = (int)(mass*1000);
       //     System.out.println(massloc);
            int massLocation = massloc -startRange;
            int cs = spectra.chargeState;
            if(indexedSpectraListTable[cs]==null)
            {
                indexedSpectraListTable[cs] = new ArrayList<>();
                massIndex[cs] = new int[endRange - startRange+100];

            }
            indexedSpectraListTable[cs].add(spectra);
            massIndex[cs][massLocation]++;
           // System.out.print("\r");
        }
        fillIndex();
        sortLibrarySpectra();

        return true;
    }


    public void getLibraryRanges(String [] libArray, String path) throws IOException {
        BufferedReader br;
        String line;
        for(String libFiles: libArray)
        {
            if(libFiles.endsWith("ms2")) {
                String hString = getRangeString(path+File.separatorChar+libFiles);
                String [] arr = hString.split("\t");
                int start  = Integer.parseInt(arr[2]);
                int end = Integer.parseInt(arr[3]);
                if(start<0) continue;
                for(int i = start; i<end; i++)
                {
                    if(fileIndex[i]==null)
                    {
                        fileIndex[i] = new ArrayList<>();
                    }
                    fileIndex[i].add(libFiles);
                }


               /* List<String> startList= massRangeFileMap.get(start);
                List<String> endList= massRangeFileMap.get(start);
                if(startList!=null)
                {
                    startList.add(hString);
                }
                if(endList!=null)
                {
                    endList.add(hString);
                }
                if(startList== null && endList == null)
                {
                    Range<Integer> range = Range.closed(start,end);
                    List<String> list = new ArrayList<>();
                    list.add()
                    massRangeFileMap.put(range,list);
                }*/
                //massRangeFileMap.put(range,libFiles);
            }
        }
    }

    public static String getRangeString(String ms2file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(ms2file));
        String line = br.readLine();
        br.close();
        return line;
    }

    /*public RangeMap<Integer, String> getMassRangeFileMap() {
        return massRangeFileMap;
    }*/

    public edu.scripps.pms.mspid.SearchParams getSearchParams()
    {
        return params;
    }

    public int getStartRange() {
        return startRange;
    }

    public int getEndRange() {
        return endRange;
    }

    public List<String> getFilesInRange()
    {

        Set<String> resultSet = new HashSet<>();
        for(int i=startRange; i<endRange; i++)
        {
            if(fileIndex[i]!=null)
            {
                resultSet.addAll(fileIndex[i]);
            }
        }

        return new ArrayList<>(resultSet);
    }

    public void setLibraryIndexer(LibraryIndexer libraryIndexer)
    {
        this.libraryIndexer = libraryIndexer;
    }


    public synchronized void close() throws IOException {
        if(bw != null) bw.close();
        bw = null;
    }

    public synchronized void write(String output) throws IOException {
        if(bw==null) bw = new BufferedWriter(new FileWriter(outputPath));
        bw.write(output);
        bw.newLine();
    }

    public String getOutputPath() {
        return outputPath;
    }

    public void setOutputPath(String outputPath) {
        this.outputPath = outputPath;
    }
}
