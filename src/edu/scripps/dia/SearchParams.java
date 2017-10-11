package edu.scripps.dia;

/**
 * Created by yateslab on 8/7/17.
 */
/**
 * @file SearchParams.java
 * This is the source file for edu.scripps.pms.mspid.SearchParams
 * @author Tao Xu
 * @date $Date: 2014/07/08 23:52:13 $
 */

import edu.scripps.pms.mspid.*;
import edu.scripps.pms.util.enzyme.Protease;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;



public class SearchParams extends edu.scripps.pms.mspid.SearchParams {
    private String databaseName;
    private int primaryScoreType = 0;
    private int secondaryScoreType = 1;
    private int minMatch;
    private String precursorIsotope;
    private String fragmentIsotope;
    private int numIsotopicPeaks;
    private int preprocess;   // preprocess mode
    private double highPrecursorTolerance;
    private double lowPrecursorTolerance;
    private double fragmentTolerance;
    private double precursorTolerance;
    private float precursorToleranceFloat;
    private double maxPrecursorMass = 10000;
    private double minPrecursorMass = 600;
    private int minimumPeptideLength = 5;
    private int maxNumSpectra = 5000;
    private int minNumSpectra = 25;
    private int maxAlter;
    private int maxPrecursorCharge = 10000;
    private int minPrecursorCharge = 0;
    private boolean isDatabaseIndexed;
    private int peakRankThreshold;
    private int candidatePeptideThreshold;
    private int numOutput = 5;
    private int locusType = 0;
    private int displayMod = 0;
    private String atomicEnrichment = "0";
    //    private int incompletelabeling = 0;
    private int leftshift1 = 1000000;
    private int leftshift2 = 1000000;
    private int leftshift3 = 1000000;
    //    private Modifications mods = new Modifications();
    private boolean isDeCharged = false;
    private boolean isEtdSearch = false;
    private int multistageActivationMod = 0;
    private String paramFilePath;
    private String paramFileName;
    private boolean chargeDisambiguation = false;
    private Protease p;
    private int enzymeSpecificity;
    private int maxInternalMisCleavage = -1; // -1 for unlimited
    private ArrayList<Modification> staticMods = new ArrayList<Modification>();
    private ArrayList<DiffMod> diffMods = new ArrayList<DiffMod>();
    private ArrayList<TerminalModification> ntermDiffMods = new ArrayList<TerminalModification>();
    private ArrayList<TerminalModification> ctermDiffMods = new ArrayList<TerminalModification>();

    private double ntermStaticMod;
    private double ctermStaticMod;
    private ArrayList<Modifications> allModifications = new ArrayList<Modifications>(1000);
    private HashSet<Modifications> mods = new HashSet<Modifications>();
    private double maxMassShift = -10000;

    public SearchParams(String path, String paramFile) throws IOException, JDOMException {
        paramFilePath = path;
        paramFileName = paramFile;
        readParams();
    }

    public SearchParams(String paramFile) throws IOException, JDOMException {
        paramFileName = paramFile;
        //paramFilePath = ".";
        readParams();
    }
    public SearchParams() {

    }
    public boolean isEtdSearch() {
        return isEtdSearch;
    }
    protected Modifications getStaticTerminalModifications() {

        Modifications temp = new Modifications();
        temp.setStaticNTermMod(ntermStaticMod);
        temp.setStaticCTermMod(ctermStaticMod);
        return temp;
    }
    public void setMaxInternalMisCleavage(int msc) {
        maxInternalMisCleavage = msc;
    }

    public int getMaxInternalMisCleavage() {
        if(enzymeSpecificity == 0) return -1; // -1 for unlimited
        return maxInternalMisCleavage;
    }

    public boolean isDeCharged() {
        return isDeCharged;
    }
    public int getMultistageActivationMod() {
        return multistageActivationMod;
    }
    public int getNumNTermDiffMods() {
        return ntermDiffMods.size();
    }
    public Iterator<TerminalModification> getNTermDiffMods() {
        return ntermDiffMods.iterator();
    }
    public int getNumCTermDiffMods() {
        return ctermDiffMods.size();
    }
    public Iterator<TerminalModification> getCTermDiffMods() {
        return ctermDiffMods.iterator();
    }
    public int getNumModifications() {
        return allModifications.size();
    }
    public double getMaxMassShift() {
        if(maxMassShift < -9999 ) {

            for(Iterator<Modifications> it = getAllModifications(); it.hasNext();) {
                double shift = it.next().getMassShift();
                maxMassShift = maxMassShift > shift? maxMassShift : shift;
            }
            maxMassShift = maxMassShift < -9999? 0 : ntermStaticMod + ctermStaticMod; // in case no modification
        }
        return maxMassShift;

    }
    public Iterator<Modifications> getAllModifications() {
        return allModifications.iterator();
    }
    public double getStaticCTermMod() {
        return ctermStaticMod;
    }
    public double getStaticNTermMod() {
        return ntermStaticMod;
    }

    public int getNumOutput() {
        return numOutput;
    }
    private void addAllModifications() {
        int numDiffMods = diffMods.size();
        Modifications m = getStaticTerminalModifications();
        addTerminalDiffMods(m); // with terminal diff mods only
        for(int diffModIndex = 0; diffModIndex < numDiffMods; diffModIndex++) {
            addDiffMods(diffModIndex, maxAlter, m);
        }
        //for(Iterator<Modifications> it = allModifications.iterator(); it.hasNext();) {
        //    System.out.println(it.next().getInfo());
        //}
    }
    private void addModifications(Modifications m) {
        //System.out.print("begin in addModifications, massShift: " + m.getDiffModsShift() + "\tnumModificationsAddeSoFar: " + allModifications.size() + "\t");

        if(m != null && m.getDiffModsShift() != 0 && !allModifications.contains(m)) {
            if(allModifications.contains(m) ||  m.getDiffModsShift() == 0) {
//System.out.println ("contains " + m.getNTermMassShift() + " " + m.getCTermMassShift());
            }
            allModifications.add(m.copy());
        } else {
            //System.out.println(m.getInfo() + "Not added!!!");
        }
    }

    // add c and n terminal modifications and then add m to allModifications
    private void addTerminalDiffMods(Modifications mToAdd) {
        Modifications m = mToAdd.copy();
        int numNtermMods = ntermDiffMods.size();
        int numCtermMods = ctermDiffMods.size();
        addModifications(m); // no n no c
        if(numNtermMods == 0) {
            if(numCtermMods == 0) {
                return;
            } else {
                for(int i = 0; i < numCtermMods; i++) {
                    Modifications newm = m.copy();
                    newm.setCTerm(ctermDiffMods.get(i));
                    addModifications(newm);
                }
            }
        } else {
            if(numCtermMods == 0) {
                for(int i = 0; i < numNtermMods; i++) {
                    Modifications newm = m.copy();
                    newm.setNTerm(ntermDiffMods.get(i));
                    addModifications(newm);
                }
            } else {
                for(int i = 0; i < numNtermMods; i++) {
                    Modifications newm = m.copy();
                    newm.setNTerm(ntermDiffMods.get(i)); // n no c
                    addModifications(newm);
                    for(int j = 0; j < numCtermMods; j++) {
                        Modifications nc = newm.copy();
                        nc.setCTerm(ctermDiffMods.get(j));
//System.out.println("in addterminal: " + nc.getNTermMassShift() + "cterm: " + nc.getCTermMassShift());
                        addModifications(nc);       // n + c
                    }
                }
                for(int i = 0; i < numCtermMods; i++) {
                    Modifications newm = m.copy();
                    newm.setCTerm(ctermDiffMods.get(i));
                    addModifications(newm);  // no n + c
                }
            }
        }
        //System.out.println("end in addTerminalDiffMods");
    }
    private void addDiffMods(int index, int numToAdd, Modifications temp) {
        if(numToAdd < 1 || index >= diffMods.size() || temp.getNumDiffMods() >= maxAlter) {
            return;
        }
        //addModifications(m);
        Modifications newm = temp.copy();
        //System.out.println("in addDiffMods, indiex: " + index + "\tnumToAdd: " + numToAdd + "\tnumModsInTemp: " + temp.getInfo());
        //addTerminalDiffMods(m); // this will add  m and its siblings ith terminal


        addDiffMods(index+1, numToAdd, newm); // without add this one
        for(int i = 1; i <= numToAdd; i++) {
            newm = newm.copy();
            newm.addDiffMod(diffMods.get(index)); // add one
            addTerminalDiffMods(newm);
            addDiffMods(index+1, numToAdd-i, newm); // after add one

        }
    }
    private void readParams() throws IOException, JDOMException {
        String file = paramFilePath == null? paramFileName : paramFilePath + "/" + paramFileName;
        //File paraFile = new File(paramFilePath + "/" + paramFileName);
        File paraFile = new File(file);
        FileInputStream paramInput = new FileInputStream(paraFile);
        Document doc = new SAXBuilder().build(paramInput);

        Element root = doc.getRootElement();
        readDatabase(root.getChild("database"));
        readSearchMode(root.getChild("search_mode"));
        readIsotopes(root.getChild("isotopes"));
        readTolerance(root.getChild("tolerance"));
        readPrecursorMassLimits(root.getChild("precursor_mass_limits"));
        readPrecursorChargeLimits(root.getChild("precursor_charge_limits"));
        readPeptideLengthLimits(root.getChild("peptide_length_limits"));

        Element numpeaklimit = root.getChild("num_peak_limits");
        numpeaklimit = numpeaklimit == null? root.getChild("num_spectra_limits") : numpeaklimit;
        if(numpeaklimit != null) {
            readNumSpectraLimits(numpeaklimit);
        }

        String maxnumdiffmods = root.getChildTextTrim("max_num_diffmod");
        maxnumdiffmods = maxnumdiffmods == null? root.getChildTextTrim("max_alter") :  maxnumdiffmods;

        maxAlter = (maxnumdiffmods == null || "".equals(maxnumdiffmods))? 0 : Integer.parseInt(maxnumdiffmods);

        readModifications(root.getChild("modifications"));
        readEnzymeInfo(root.getChild("enzyme_info"));

        paramInput.close();
    }
    private void readEnzymeInfo(Element e) {
        enzymeSpecificity = Integer.parseInt(e.getChildTextTrim("specificity"));
        String maxmiscleavage =  e.getChildTextTrim("max_num_internal_mis_cleavage");

        System.out.println("Maximum number of internal mis : " + maxmiscleavage);
        if(maxmiscleavage != null && !"".equals(maxmiscleavage.trim())) {
            maxInternalMisCleavage = Integer.parseInt(maxmiscleavage);
        }
        System.out.println("maxInternalMisCleavage value is now: " + maxInternalMisCleavage);
        String enzymeName = e.getChildTextTrim("name");
        p = new Protease(enzymeName);
        p.setType(Boolean.parseBoolean(e.getChildTextTrim("type")));
        Element re = e.getChild("residues");
        List residues = re.getChildren();
        for(Iterator i = residues.iterator(); i.hasNext();) {
            Element r = (Element)i.next();
            String s = r.getTextTrim();
            p.addCleavageSite(s.charAt(0));
        }

    }
    private void readDatabase(Element e) {
        databaseName = e.getChildTextTrim("database_name");
        String s = e.getChildTextTrim("is_indexed");
        if(s != null) {
            char c = s.charAt(0);
            if(c == 'T' || c == 't' || c == 'Y' || c == 'y') {
                isDatabaseIndexed = true;
            } else if(c == 'N' || c == 'n' || c == 'F' || c == 'f') {
                isDatabaseIndexed = false;
            } else {
                throw new InvalidArgumentException(
                        "is_indexed unknown. Should be true or false, or yes or no.");
            }
        } else {
            isDatabaseIndexed = false; // by default, database is not indexed
        }
    }
    private void readTerminalMods(Element me, String terminal) {
        boolean isCterm = terminal.startsWith("c");
        Element e = me.getChild(terminal);
        if(e != null) {
            Element staticMod = e.getChild("static_mod");
            if(staticMod != null) {
//System.out.println("in readTerminalMods, terminal: " + terminal + "\t" + staticMod.getChildTextTrim("mass_shift"));
                //double massShift = Double.parseDouble(staticMod.getChildTextTrim("mass_shift"));
                double massShift = 0;
                String temp = staticMod.getChildTextTrim("mass_shift");
                if(temp != null && !"".equals(temp.trim())) {
                    massShift = Double.parseDouble(temp);
                }
                if(isCterm) {
                    ctermStaticMod = massShift;
                } else {
                    ntermStaticMod = massShift;
                }
                //    char symbol = e.getChildTextTrim("symbol").charAt(0);
            }
//System.out.println("in readTerminalMods, ctermStaticMod: " + ctermStaticMod + "\tntermStaticMod: " + ntermStaticMod);
            Element diffMods = e.getChild("diff_mods");
            if(diffMods != null) {
                List <Element> diffModList = diffMods.getChildren();
                for(Element de : diffModList) {
                    //char residue = e.getChildTextTrim("residue").charAt(0);
                    char symbol = de.getChildTextTrim("symbol").charAt(0);
                    String massShift = de.getChildTextTrim("mass_shift");
                    //System.out.println("static modification: " + staticModification);
                    if (massShift != null && !"".equals(massShift.trim())) {
                        //m.setStaticModification(Double.parseDouble(staticModification));
                        double mass = Double.parseDouble(massShift);
                        if(mass != 0) {
                            TerminalModification d = new TerminalModification(symbol, mass);
                            if(isCterm) {
                                ctermDiffMods.add(d);
                            } else {
                                ntermDiffMods.add(d);
                            }
                            if(displayMod == 0) {
                                d.setModInfo("(" + d.getMassShift() + ")");
                            } else {
                                d.setModInfo("" + symbol);
                            }
                        }
                        //System.out.println(terminal + " diff mod: " + massShift);
                    }
                }

            }
        }
    }
    private void readModifications(Element me) {
        String moddisplay = me.getChildTextTrim("display_mod");
        if(moddisplay != null) {
            displayMod = Integer.parseInt(moddisplay);
        }

        //System.out.println("In readModifications, moddisplay: " + moddisplay + " and display_mod: " + displayMod);
        readTerminalMods(me, "c_term");
        readTerminalMods(me, "n_term");


        // get static mods
        List <Element> staticModList = me.getChild("static_mods").getChildren();
        for(Element e : staticModList) {
            char residue = e.getChildTextTrim("residue").charAt(0);
            String massShift = e.getChildTextTrim("mass_shift");
            //System.out.println("static modification: " + staticModification);
            if (massShift != null) {
                double mass = Double.parseDouble(massShift);
                if(mass != 0) {
                    staticMods.add(new Modification(residue, mass));
                }
                //System.out.println("static modification: on " + residue + "\t" +  massShift);
            }
        }
        // get diff_mods
        List <Element> diffModList = me.getChild("diff_mods").getChildren();
        for(Element e : diffModList) {
            //char residue = e.getChildTextTrim("residue").charAt(0);
            char symbol = e.getChildTextTrim("symbol").charAt(0);
            String massShift = e.getChildTextTrim("mass_shift");

            //System.out.println("static modification: " + staticModification);
            if (massShift != null) {
                //m.setStaticModification(Double.parseDouble(staticModification));
                double mass = Double.parseDouble(massShift);
                if(mass != 0) {
                    String neutralLoss = e.getChildTextTrim("neutral_loss");
                    double nl = 0;
                    if(neutralLoss != null) {
                        nl = Double.parseDouble(neutralLoss);
                    }
                    StringBuffer sb = new StringBuffer(200);
                    sb.append("Diff mod: " + mass + "\tsymbol: " + symbol + "\tresidue:");
                    DiffMod d = new DiffMod(mass, nl, symbol);
                    if(displayMod == 0) {
                        d.setModInfo("(" + d.getMassShift() + ")");
                    } else {
                        d.setModInfo(""+symbol);
                    }
                    List<Element> residues = e.getChild("residues").getChildren();
                    for(Element r : residues) {
                        byte c = (byte)r.getTextTrim().charAt(0);
                        sb.append("\t" + (char)c);
                        d.setModifiable(c, true);
                    }
                    diffMods.add(d);
                    //System.out.println(sb);
                    //System.out.println("info: " + d.getModInfo());
                }
            }
        }
        addAllModifications();

        if(maxAlter > 0) {
            System.out.println("Number of Internal Diff Modification types: " + diffMods.size());
            System.out.println("Total Number of Differential Modifications (mass shift) considered: " + getNumModifications());
        }
    }
    private void readNumSpectraLimits(Element ne) {
        if(ne == null) return;
        maxNumSpectra =  Integer.parseInt(ne.getChildTextTrim("maximum"));
        minNumSpectra =  Integer.parseInt(ne.getChildTextTrim("minimum"));
    }
    private void readPrecursorMassLimits(Element pe) {
        if(pe == null) return;
//        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
//        minPrecursorMass = Double.parseDouble(pe.getChildTextTrim("minimum"));
        minPrecursorMass = Double.parseDouble(pe.getChildTextTrim("minimum"));
    }
    private void readPrecursorChargeLimits(Element pe) {
//        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
        if(pe == null) {
            return;
        }
        maxPrecursorCharge = Integer.parseInt(pe.getChildTextTrim("maximum"));
        minPrecursorCharge = Integer.parseInt(pe.getChildTextTrim("minimum"));
    }

    private void readPeptideLengthLimits(Element pe) {
        if(pe == null) return;
        minimumPeptideLength = Integer.parseInt(pe.getChildTextTrim("minimum"));
    }

    private void readTolerance(Element te) {
//        precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        //precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        highPrecursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor_high"));
        lowPrecursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor_low"));
//        fragmentTolerance = Double.parseDouble(te.getChildTextTrim("fragment"));
        precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        this.precursorToleranceFloat = (float)precursorTolerance/1000;




        fragmentTolerance = Double.parseDouble(te.getChildTextTrim("fragment"));

    }
    private void readIsotopes(Element ie)
            throws InvalidArgumentException {
        String mono = MassSpecConstants.MONOISOTOPE;
        String avg = MassSpecConstants.AVGISOTOPE;
        String isotope =  ie.getChildTextTrim("precursor");
        numIsotopicPeaks = Integer.parseInt(ie.getChildTextTrim("num_peaks"));
        // set the precursor isotope type
        if (mono.equals(isotope) || avg.equals(isotope)) {
            precursorIsotope = isotope;
        } else {
            throw new InvalidArgumentException(
                    "Unknown precursor isotope type. Should be mono or avg.");

        }

        // set fgragment isotope type
        isotope = ie.getChildTextTrim("fragment");
        if (mono.equals(isotope) || avg.equals(isotope)) {
            fragmentIsotope =  isotope;
        } else {
            throw new InvalidArgumentException(
                    "Unknown fragment isotope type. Should be mono or avg.");
        }
    }
    private void readSearchMode(Element me) {
        String primaryscoretype = me.getChildTextTrim("primary_score_type");
        String secondaryscoretype = me.getChildTextTrim("secondary_score_type");
        primaryScoreType = primaryscoretype == null? 0 : Integer.parseInt(primaryscoretype);
        if(secondaryscoretype != null) {
            secondaryScoreType = Integer.parseInt(secondaryscoretype);
        } else {
            secondaryScoreType = primaryScoreType == 0? 1 : 0;
        }
        String peakrankt = me.getChildTextTrim("peak_rank_threshold");
        peakRankThreshold =  peakrankt == null? 200: Integer.parseInt(peakrankt);

        String candidatet = me.getChildTextTrim("candidate_peptide_threshold");
        candidatePeptideThreshold = candidatet == null? 500: Integer.parseInt(candidatet);

        String numoutput = me.getChildTextTrim("num_output");
        numOutput = numoutput == null? 5 : Integer.parseInt(numoutput);
        minMatch =  Integer.parseInt(me.getChildTextTrim("min_match"));

        String preproc = me.getChildTextTrim("preprocess");
        preprocess  =  preproc == null? 1 : Integer.parseInt(preproc);
        String fragMethod = me.getChildTextTrim("fragmentation_method");
        isEtdSearch = fragMethod != null && (fragMethod.startsWith("ETD") || fragMethod.startsWith("etd"));

        String locus = me.getChildTextTrim("locus_type");
        locusType = locus == null? 0 : Integer.parseInt(locus);

        atomicEnrichment = me.getChildTextTrim("atomic_enrichement");
        String leftshift = atomicEnrichment;
        if(leftshift != null) {
            //System.out.println("Atomic enrichement: " + leftshift);

            double atomicEnrichment = Double.parseDouble(leftshift);
            int ae = 0;
            if(atomicEnrichment != 0) {
                ae = (int)(atomicEnrichment);
            }

            switch(ae) {
                case 90 : leftshift1 = 2; leftshift2 = 4; leftshift3 = 8; break;
                case 91 : leftshift1 = 2; leftshift2 = 5; leftshift3 = 9; break;
                case 92 : leftshift1 = 2; leftshift2 = 6; leftshift3 = 10; break;
                case 93 : leftshift1 = 2; leftshift2 = 6; leftshift3 = 11; break;
                case 94 : leftshift1 = 3; leftshift2 = 8; leftshift3 = 14; break;
                case 95 : leftshift1 = 3; leftshift2 = 9; leftshift3 = 16; break;
                case 96 : leftshift1 = 4; leftshift2 = 12; leftshift3 = 21; break;
                case 97 : leftshift1 = 5; leftshift2 = 16; leftshift3 = 29; break;
                case 98 : leftshift1 = 7; leftshift2 = 24; leftshift3 = 36; break;
                case 99 : leftshift1 = 10; leftshift2 = 36; leftshift3 = 1000; break;

            }
        }
        //System.out.println("leftshift1: " + leftshift1 + "\tleftshift2: " + leftshift2 + "\tleftshift3: " + leftshift3);
        String chargedisamb = me.getChildTextTrim("charge_disambiguation");
        if(chargedisamb != null && Integer.parseInt(chargedisamb) != 0) {
            chargeDisambiguation = true;
        }

        String decharged = me.getChildTextTrim("is_decharged");
        int dechargestate = decharged == null? 0 : Integer.parseInt(decharged);
        isDeCharged = dechargestate == 1? true : false;

        String multistageActivation = me.getChildTextTrim("multistage_activation_mode");
        multistageActivationMod = multistageActivation == null? 0 : Integer.parseInt(multistageActivation);
    }
    public boolean getChargeDisambiguation() {

        return chargeDisambiguation;
    }
    public int getLeftShift1() {
        return leftshift1;
    }
    public int getLeftShift2() {
        return leftshift2;
    }
    public int getLeftShift3() {
        return leftshift3;
    }
    public void setDbName(String dbName) {
        databaseName = dbName;
    }
    public void setPrimaryScoreType(int primaryScoreType) {
        this.primaryScoreType = primaryScoreType;
    }

    public void setsecondaryScoreType(int secondaryScoreType) {
        this.secondaryScoreType = secondaryScoreType;
    }
    public void setMinMatch(int minMatch) {
        this.minMatch = minMatch;
    }

    public void setPrecursorIsotope(String precursorIsotope)
            throws InvalidArgumentException {

        if (MassSpecConstants.MONOISOTOPE.equals(precursorIsotope) ||
                MassSpecConstants.AVGISOTOPE.equals(precursorIsotope)) {
            this.precursorIsotope = precursorIsotope;
        } else {
            throw new InvalidArgumentException(
                    "Unknow precursor isotope type. Should be mono or avg.");
        }
    }
    public void setFragmentIsotope(String fragmentIsotope)
            throws InvalidArgumentException {

        if (MassSpecConstants.MONOISOTOPE.equals(fragmentIsotope) ||
                MassSpecConstants.AVGISOTOPE.equals(fragmentIsotope)) {
            this.fragmentIsotope = fragmentIsotope;
        } else {
            throw new InvalidArgumentException(
                    "Unknow fragment isotope type. Should be mono or avg.");
        }
    }

    public void setPreprocess(int mode) {
        preprocess = mode;
    }
    public void setFragmentTolerance(double tolerance) {
        fragmentTolerance = tolerance;
    }
    public void setMaxPrecursorMass(double mass) {
        maxPrecursorMass = mass;
    }
    public void setMinPrecursorMass(double mass) {
        minPrecursorMass = mass;
    }
    public void setMaxPrecursorCharge(int charge) {
        maxPrecursorCharge = charge;
    }
    public void setMinPrecursorCharge(int charge) {
        minPrecursorCharge = charge;
    }
    public void setMinimumPeptideLength(int len) {
        minimumPeptideLength = len;
    }
    public void setMaxNumSpectra(int numSpectra) {
        maxNumSpectra = numSpectra;
    }
    public void setMinNumSpectra(int numSpectra) {
        minNumSpectra = numSpectra;
    }
    public void setMaxAlter(int maxAlter) {
        this.maxAlter = maxAlter;
    }
    public int getLocusType() {
        return locusType;
    }
    public int getPeakRankThreshold() {
        return peakRankThreshold;
    }
    public int getCandidatePeptideThreshold() {
        return candidatePeptideThreshold;
    }
    public Iterator<Modification> getStaticMods() {
        return staticMods.iterator();
    }
    public Iterator<DiffMod> getDiffMods() {
        return diffMods.iterator();
    }
    public int getNumDiffMods() {
        return diffMods.size();
    }
    public String getDbName() {
        return databaseName;
    }
    public boolean isDatabaseIndexed() {
        return isDatabaseIndexed;
    }

    public int getPrimaryScoreType() {
        return primaryScoreType;
    }

    public int getSecondaryScoreType() {
        return secondaryScoreType;
    }

    public int getMinMatch() {
        return minMatch;
    }
    public int getNumIsotopicPeaks() {
        return numIsotopicPeaks;
    }
    public String getPrecursorIsotope() {
        return precursorIsotope;
    }

    public String getFragmentIsotope() {
        return fragmentIsotope;
    }
    public int getPreprocess() {
        return preprocess;
    }

    public double getHighPrecursorTolerance() {
        return highPrecursorTolerance;
    }
    public double getLowPrecursorTolerance() {
        return lowPrecursorTolerance;
    }
    public double getPrecursorTolerance() {
        return precursorTolerance;
    }
    public double getFragmentTolerance() {
        return fragmentTolerance;
    }

    public int getMinimumPeptideLength() {
        return minimumPeptideLength;
    }
    public double getMaxPrecursorMass() {
        return maxPrecursorMass;
    }
    public double getMinPrecursorMass() {
        return minPrecursorMass;
    }
    public double getMaxPrecursorCharge() {
        return maxPrecursorCharge;
    }
    public double getMinPrecursorCharge() {
        return minPrecursorCharge;
    }

    public int getMaxNumSpectra() {
        return maxNumSpectra;
    }
    public int getMinNumSpectra() {
        return minNumSpectra;
    }
    public int getMaxAlter() {
        return maxAlter;
    }
    public Protease getProtease() {
        return p;
    }
    public int getEnzymeSpecificity() {
        return enzymeSpecificity;
    }
    public String getEnzymeName() {
        return p.getName();
    }
    public Protease getEnzyme() {
        return p;
    }
    public double getStaticTerminalMods() {
        //return mods.getStaticTerminalMods();
        return ctermStaticMod + ntermStaticMod;
    }
    public String getDatabaseName() {
        return databaseName;
    }
    public int getDisplayMod () {
        return displayMod ;
    }
    public String getAtomicEnrichment () {
        return atomicEnrichment;
    }
    public void setDatabaseName(String dbname) {
        databaseName = dbname;
    }
    public String output() {
        StringBuffer sb = new StringBuffer(10000);
        sb.append("<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n");
        sb.append("<!--  Parameters for ProLuCID Seach  -->\n");
        sb.append("<parameters>\n");
        sb.append("\t<database>\n");
        sb.append("\t\t<database_name>" + getDatabaseName() + "</database_name>\n");
        sb.append("\t</database>\n");

        sb.append("</parameters>\n");
        return sb.toString();
    }
    public static void main(String args[]) throws Exception {
        /*
        SearchParams sp = new SearchParams(args[0], args[1]);
        System.out.println("Database is: " + sp.getDbName());
        System.out.println("maxAlter is: " + sp.getMaxAlter());
        System.out.println("highPrecursorTolerance is: " + sp.getHighPrecursorTolerance());
        System.out.println("Fragment isotope is: " + sp.getFragmentIsotope());
        Iterator <Modification> modifications = sp.getStaticMods();
        while(modifications.hasNext()) {
            Modification m = modifications.next();
            System.out.println("Residue " + (char)m.getResidue() + "\t" + m.getMassShift() + "\n");
        }

        String a = "java";
        String b = new StringBuffer(a).toString();
        a = b;
        if(a.equals(b)) {
            System.out.println("a equal b");

        }
        if(a == b) {

            System.out.println("a == b");
        }else {

            System.out.println("a != b");
        }
        */
        edu.scripps.pms.mspid.SearchParams sp = new edu.scripps.pms.mspid.SearchParams("/home/rpark/search.xml");

        System.out.println(sp.output());
        System.out.println(sp.getDbName());
    }

    public float getPrecursorToleranceFloat() {
        return precursorToleranceFloat;
    }

    public void setPrecursorToleranceFloat(float precursorToleranceFloat) {
        this.precursorToleranceFloat = precursorToleranceFloat;
    }
}




