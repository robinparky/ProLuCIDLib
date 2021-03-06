package edu.scripps.dia;

import edu.scripps.dia.util.MzxmlPeakList;
import edu.scripps.dia.util.MzxmlSpectrumReader;
import edu.scripps.dia.util.SpectrumReader;
import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;
import gnu.trove.TIntArrayList;
import gnu.trove.TIntHashSet;
import gnu.trove.TIntObjectHashMap;
import org.jdom.JDOMException;

import java.io.*;
import java.sql.SQLException;
import java.util.*;

/**
 * Created by yateslab on 8/8/17.
 */
public class ConsensusLibrarySearchEngine {
    private static int NUMFINALRESULT = 5;

    private int tempMaxCS = 20;
    private List<PeakList> peakLists;
    private List<PeakList> libraryPeakLists = new ArrayList<>();
    private List<ProcessedPeakList>[] libraryPeakListTable = new List[tempMaxCS];
    private int [][] massIndex = new int[tempMaxCS][];
    private SearchParams params;
    private String ms2FilePath;
    private LibrarySQLiteManager library;
    private MassCalculator mc;
    private int startRange;
    private int endRange;
    private boolean useSQLite = false;
    private static DistributionCalculator dc;

    public static void main(String[] args) throws Exception {
        String ms2Path = args[0];
        String paramsPath = args[1];
        String libraryPath = args[2];
        String output = args[3];
        ConsensusLibrarySearchEngine clse = new ConsensusLibrarySearchEngine(ms2Path, paramsPath);
        clse.calcRange();
        clse.readMS2TargetFile(libraryPath);
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


    public ConsensusLibrarySearchEngine(String ms2Path, String searchParamsPath) throws IOException, JDOMException, SQLException {
        params = new SearchParams(searchParamsPath);
        ms2FilePath = ms2Path;
        mc = new MassCalculator(params);
        if(dc == null) {
            dc = new DistributionCalculator();
        }
    }

    public ConsensusLibrarySearchEngine(String ms2Path, String searchParamsPath, String dbPath) throws IOException, JDOMException, SQLException {
        params = new SearchParams(searchParamsPath);
        ms2FilePath = ms2Path;
        library = new LibrarySQLiteManager(dbPath);
        mc = new MassCalculator(params);
        useSQLite = false;
    }

    public ConsensusLibrarySearchEngine(String ms2Path, SearchParams params) {
        this.params = params;
        ms2FilePath = ms2Path;
    }

    private void init() throws Exception {
        peakLists = getSpectra(ms2FilePath);
    }

    public SearchResult search(PeakList peakList, Zline zline) throws SQLException {
        int numPeaks = peakList.numPeaks();
        PriorityQueue<ScoredPeptideHit> sphQueue = new PriorityQueue<>(new Comparator<ScoredPeptideHit>() {
            @Override
            public int compare(ScoredPeptideHit scoredPeptideHit, ScoredPeptideHit t1) {
                return Double.compare(scoredPeptideHit.getPScore(),t1.getPScore());
            }
        });
        TIntHashSet set = new TIntHashSet();
        double worseScore = Double.MAX_VALUE;
        SearchResult searchResult = null;
        System.out.println(peakList.getLoscan());
        if (numPeaks > params.getMinNumSpectra() && numPeaks < params.getMaxNumSpectra()) {
            float prcMass = (float) zline.getM2z();
            ScoredPeptideHit sph;
            if (prcMass > params.getMinPrecursorMass() && prcMass < params.getMaxPrecursorMass()) {
                TIntArrayList highLimits = new TIntArrayList();
                TIntArrayList lowLimits = new TIntArrayList();
                MassRangeFinder.findRange(prcMass, params, lowLimits, highLimits);
                ProcessedPeakList ppl = new ProcessedPeakList(peakList, zline, params, mc);
                ppl.preprocess(params.getPreprocess());
                searchResult = new SearchResult(ppl);
                set.clear();
                for (int i = 0; i < highLimits.size(); i++) {
                    int high = highLimits.get(i);
                    int low = lowLimits.get(i);
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
                        // sph = ppl.prelimScoreCorrelation(pl);
                        System.out.println(pl.getLoscan()+"\t"+sph.getPScore());
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
            }
        }
        System.out.println();
        //searchResult.setFinalResultsForPrelimScores();
        List<ScoredPeptideHit> sphList = new ArrayList<>(sphQueue);
        int size = sphList.size()<NUMFINALRESULT?sphList.size():NUMFINALRESULT;
        searchResult.setPrelimScoreHits(sphList,size);
        return searchResult;
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


    private static List<PeakList> getSpectra(String file) throws IOException, JDOMException, Exception {
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

    public void readMS2TargetFile(String ms2FilePath) throws IOException {
        SpectrumReader reader = new SpectrumReader(ms2FilePath,"ms2");
        List<PeakList> plList = reader.getSpectraList();
        reader.closeDataFile();
        for(PeakList peakList:plList )
        {
           int z=  peakList.getFirstChargeState();
           peakList.setFilename(ms2FilePath);
           double mass = peakList.getZlines().next().getM2z();
           int massLocation = (int) (mass*1000) -startRange;
           if(libraryPeakListTable[z]==null)
           {
               libraryPeakListTable[z] = new ArrayList<>();
               massIndex[z] = new int[endRange - startRange];
           }
           ProcessedPeakList ppl = new ProcessedPeakList(peakList,peakList.getZlines().next(),params,mc);
           libraryPeakListTable[z].add(ppl);
           massIndex[z][massLocation]++;
        }
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
        int startrange=-1;
        int endrange=-1;
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

            System.out.println("" + ms2FilePath);
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
        this.endRange = endrange;
        br.close();
    }

}
