/**
 * @file ModifiedPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.ModifiedPeptideHit
 * @author Tao Xu
 * @date $Date: 2010/03/05 23:34:55 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;

import java.util.ArrayList;
import java.util.Iterator;

public class ModifiedPeptideHit extends PeptideHit {
     
    protected Modifications mods;
    protected String extendedSeq = null;
    protected String exactSeq = null; // same as extendedSeq except leading and tailing residue
    protected DiffMod [] diffMods;
    //protected ArrayList<ModifiedPeptideHit> candidates;
    protected ArrayList<DiffModSite> sites = new ArrayList<DiffModSite>();
    protected int numPotentialSites = -1;
    protected double internalDiffModMassAdded = 0;
    protected int numInternalDiffModsAdded = 0;
    public static final boolean isModified = true;



    public ModifiedPeptideHit(Fasta parent, int start, int end, Modifications m) {
        super(parent, start, end);
        this.mods = m;
        diffMods = new DiffMod[end-start+1];
    }  
    public ModifiedPeptideHit(Peptide p, Modifications m) {
        super(p);
        mods = m;
        diffMods = new DiffMod[p.getEnd() - p.getStart() + 1];

    }
    public ModifiedPeptideHit(ModifiedPeptideHit mph) {
        super(mph.peptide);
        mods = mph.mods;
        diffMods = mph.diffMods;

    }
    public ModifiedPeptideHit(Peptide p) {
        super(p);
        mods = new Modifications();
        diffMods = new DiffMod[p.getEnd() - p.getStart() + 1];

    }
    public boolean isModified() {
        return isModified;
    }
    public void setDiffMod(int index, DiffMod diffmod) {
        diffMods[index] = diffmod;
        if(diffmod != null) {
            internalDiffModMassAdded += diffmod.getMassShift();
            numInternalDiffModsAdded++;
        }
    }
    public Modifications getModifications()
    {
        return mods;
    }

    public boolean isGoodCandidate() {
        return internalDiffModMassAdded == mods.getInternalDiffModsMass() &&
              numInternalDiffModsAdded == mods.getNumDiffMods();
    }
    public void setDiffMods(DiffMod [] dms) {
        this.diffMods = dms;
    }
    public DiffMod [] getDiffMods() { return diffMods;}
    protected ArrayList<ModifiedPeptideHit> generateCandidates() {
        ArrayList<ModifiedPeptideHit> candidates = new ArrayList<ModifiedPeptideHit>();
        calcModificationSites();
//System.out.println("numSites: " + numPotentialSites + "\t numDiffMods: " + mods.getNumDiffMods());
        if(getNumPotentialSites() < mods.getNumDiffMods()) {
            return candidates;
        }
        candidates.add(this.copy()); 
        int numSites = sites.size(); 

//System.out.print("numPotentialSitesInPeptide: " + numPotentialSites +  "\t numDiffMods: " + mods.getNumDiffMods());
        // may have more efficient way if numSites == mods.getNumDiffMods()
        for(int i = 0; i < numSites; i++) {
           
            DiffModSite dms = sites.get(i); 
            processDiffModSite(dms, candidates);
        }
        return getFinalCandidates(candidates);
    }
    private ArrayList<ModifiedPeptideHit> getFinalCandidates(ArrayList<ModifiedPeptideHit> candidates) {
        
        ArrayList<ModifiedPeptideHit> temp = candidates;
        candidates = new ArrayList<ModifiedPeptideHit>();
        for(Iterator<ModifiedPeptideHit> it = temp.iterator(); it.hasNext();) {
            ModifiedPeptideHit m = it.next();
            if(m.isGoodCandidate()) {
                candidates.add(m);
            }
        } 
//System.out.println("\tNumTotal Candidates: " + temp.size() + "\tNumGoodCandidates: " + candidates.size());
        return candidates;
    }
    private void processDiffModSite(DiffModSite dms, ArrayList<ModifiedPeptideHit> candidates) {
//System.out.println("numDiffModsForThisSite: " + dms.numDiffMods());
        int num = candidates.size();
        int position = dms.getPosition();
        for(int i = 0; i < num; i++) {
            ModifiedPeptideHit m = candidates.get(i);
            for(Iterator<DiffMod> it = dms.getDiffMods(); it.hasNext();) {
                DiffMod d = it.next();
                if(m.numInternalDiffModsAdded+1 <= mods.getNumDiffMods() &&
                          m.internalDiffModMassAdded+d.getMassShift() <= mods.getInternalDiffModsMass()) {
                    ModifiedPeptideHit newm = m.copy();
                    newm.setDiffMod(position, d);
                    candidates.add(newm); 
                }
            }
        }
    }   
    private ModifiedPeptideHit copy() {
        ModifiedPeptideHit newm = new ModifiedPeptideHit(peptide, mods);
        for(int i = 0; i < newm.diffMods.length; i++) {
            newm.setDiffMod(i, diffMods[i]);
        }              
        //newm.internalDiffModMassAdded = internalDiffModMassAdded;
        return newm;
    }
    protected void calcModificationSites() {
        if(numPotentialSites > -1) {
            return;
        }
        numPotentialSites = 0;
        int start = getStart();
        int end = getEnd();
        byte [] bytes = getParent().getSequenceAsBytes();
        for(int i = start; i <= end; i++) {
            if(mods.isModifiable(bytes[i])) {
                numPotentialSites++;
                sites.add(new DiffModSite(i - start, mods.residue2DiffMods(bytes[i])));
            }
        }
    }
    public int getNumPotentialSites() {
        if(numPotentialSites > -1) {
            calcModificationSites();
        }
        return numPotentialSites; 
    }    

    public Iterator<ModifiedPeptideHit> getAllModifiedPeptideHits() {
        return generateCandidates().iterator(); 
        //return candidates.iterator();
    }
    public double getTheorMass() {

        double m = ppl.getMassCalculator().getPrecursorMass(peptide.getSequence());
        m += mods.getMassShift() + mods.getNTermMassShift() + mods.getCTermMassShift();

        //return ppl.getMassCalculator().getPrecursorMass(peptide.getSequence()) +  mods.getMassShift();
        return m;
    }
    // override the function in super class
    public double getNTermStart() { // b == n term
        if(mods != null) {
            return ppl.getNTermStart() + mods.getDiffNTermMod();
        } else {
            return ppl.getNTermStart();

        }
    }    
    // override the function in super class
    public double getNTermDbinwidthStart() { // b == n term
        if(mods != null) {
            return (ppl.getNTermStart() + mods.getDiffNTermMod())*MassSpecConstants.DBINWIDTH;
        } else {
            return ppl.getNTermStart()*MassSpecConstants.DBINWIDTH;
        }
    }    
    public double getCTermStart() { // y == c term
//System.out.println("in getCTermStart in ModifiedPeptideHit, " + mods.getDiffCTermMod());
        if(mods != null) {
            return ppl.getCTermStart() + mods.getDiffCTermMod();
        } else {     
            return ppl.getCTermStart();
        }
    }    
    public double getCTermDbinwidthStart() { // y == c term
        if(mods != null) {
            return (ppl.getCTermStart() + mods.getDiffCTermMod())*MassSpecConstants.DBINWIDTH;
        } else {

            return ppl.getCTermStart()*MassSpecConstants.DBINWIDTH;
        }
    }    

    // return sequence without leading and tailing residues
    public String getExactSequence() {
        if(exactSeq != null) {
            return exactSeq;
        }
        StringBuffer sb = new StringBuffer(getLength() + 15);
        int start = peptide.getStart();
        int end = peptide.getEnd();

        double nTermMod = mods == null? 0 : mods.getDiffNTermMod();
        double cTermMod = mods == null? 0 : mods.getDiffCTermMod();
        if(nTermMod != 0) {
            sb.append("(" + nTermMod + ")"); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != 0) {
            sb.append("(" + cTermMod + ")"); 
        }

        exactSeq = sb.toString(); 
        
        return exactSeq;
    }

    public String getExtraExtendedSequence() {
        StringBuffer sb = new StringBuffer(getLength() + 25);
        // 4 more bytes for the heading and tailing aa and '.'
        int len = peptide.getLength() + 4;
        int start = peptide.getStart();
        int end = peptide.getEnd();

        switch(start){
            case 0: 
                sb.append("---");
                break;
            case 1: 
                sb.append("--");
                sb.append((char)peptide.getParent().byteAt(start-1));
                break;
            case 2: 
                sb.append('-');
                sb.append((char)peptide.getParent().byteAt(start-2));
                sb.append((char)peptide.getParent().byteAt(start-1));
                break;
            default: 
                sb.append((char)peptide.getParent().byteAt(start-3));
                sb.append((char)peptide.getParent().byteAt(start-2));
                sb.append((char)peptide.getParent().byteAt(start-1));
                break;

        }
        sb.append('.');

        TerminalModification nTermMod = mods.getNTerm();
        TerminalModification cTermMod = mods.getCTerm();
        if(nTermMod != null) {
            //sb.append("(" + nTermMod + ")"); 
            sb.append( nTermMod.getModInfo()); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != null) {
            //sb.append("(" + cTermMod + ")"); 
            sb.append(cTermMod.getModInfo()); 
        }
        sb.append('.');

        int lastIndex = peptide.getParent().getLength() - 1;

        int diff = lastIndex - end;


        switch(diff){
            case 0: 
                sb.append("---");
                break;
            case 1: 
                sb.append((char)peptide.getParent().byteAt(end+1));
                sb.append("--");
                break;
            case 2: 
                sb.append((char)peptide.getParent().byteAt(end+1));
                sb.append((char)peptide.getParent().byteAt(end+2));
                sb.append("-");
                break;
            default: 
                sb.append((char)peptide.getParent().byteAt(end+1));
                sb.append((char)peptide.getParent().byteAt(end+2));
                sb.append((char)peptide.getParent().byteAt(end+3));
                break;

        }

        return sb.toString(); 
        
        
    }
    public String getExtendedSequence() {
        if(extendedSeq != null) {
            return extendedSeq;
        }
        StringBuffer sb = new StringBuffer(getLength() + 15);
        // 4 more bytes for the heading and tailing aa and '.'
        int len = peptide.getLength() + 4;
        int start = peptide.getStart();
        int end = peptide.getEnd();
        if (start == 0) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(start-1));
        }
        sb.append('.');

        TerminalModification nTermMod = mods.getNTerm();
        TerminalModification cTermMod = mods.getCTerm();
        if(nTermMod != null) {
            //sb.append("(" + nTermMod + ")"); 
            sb.append( nTermMod.getModInfo()); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != null) {
            //sb.append("(" + cTermMod + ")"); 
            sb.append(cTermMod.getModInfo()); 
        }
        sb.append('.');

        int lastIndex = peptide.getParent().getLength() - 1;
        if (end == lastIndex) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(end+1));
        }

        extendedSeq = sb.toString(); 
        
        return extendedSeq;
    }

    public void getNLChargeNNumPeaksMatched(int z) {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        // y is used for both y and z, ie cterm ions
        // b is used for both b andc, ie nterm ions
        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function 
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        double yneutralloss = 0;
        double bneutralloss = 0;

        for(int i = start; i < end; i++) {


            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

//System.out.println((char)seq[i] + "\t" + bMass + "\t" + yMass);
            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/2f);

            // doubly charged fragment ions
            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;
 
            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            for(int j = 3; j <= z; j++) {
                int nb = (int)(0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/j);
                int ny = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/j);
                if(nb < lastTrue && nb > firstTrue) { 
                    numTheroticPeaks++;
                    if(boolMasses[nb]) numPeaksMatched++;
                }
                if(ny < lastTrue && ny > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[ny]) numPeaksMatched++;
                }
            }

            if(bneutralloss == 0 && bd != null) bneutralloss = bd.getNeutralLoss(); 
            if(yneutralloss == 0 && yd != null) yneutralloss = yd.getNeutralLoss(); 

            if(bneutralloss != 0) {
                double nlbMass = bMass - bneutralloss;
                int int2nlB = (int)(nlbMass*ppl.PRECISIONFACTOR + 0.5f);
                int int3nlB = (int) (0.5f + ppl.PRECISIONFACTOR*(nlbMass + MassSpecConstants.MASSPROTON)/2f);
                if(int3nlB < lastTrue && int3nlB > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int3nlB]) numPeaksMatched++;
    
                }
                if(int2nlB < lastTrue && int2nlB > firstTrue) {
                    //numTheroticPeaks += 2;
                    numTheroticPeaks++;
                    if(boolMasses[int2nlB]) numPeaksMatched++;
                    //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
                }
                for(int j = 3; j <= z; j++) {
                    int nlnb = (int)(0.5f + ppl.PRECISIONFACTOR*(nlbMass + MassSpecConstants.MASSPROTON)/j);
                    if(nlnb < lastTrue && nlnb > firstTrue) { 
                        numTheroticPeaks++;
                        if(boolMasses[nlnb]) numPeaksMatched++;
                    }
                }
            }
            if(yneutralloss != 0) { 
                double nlyMass = yMass - yneutralloss;
                int int2nlY = (int)(nlyMass*ppl.PRECISIONFACTOR + 0.5f);
                int int3nlY = (int)(0.5f + ppl.PRECISIONFACTOR*(nlyMass + MassSpecConstants.MASSPROTON)/2f);
                if(int3nlY < lastTrue && int3nlY > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int3nlY]) numPeaksMatched++;
                }
                if(int2nlY < lastTrue && int2nlY > firstTrue) {
                    //numTheroticPeaks += 2;
                    numTheroticPeaks++;
                    if(boolMasses[int2nlY]) numPeaksMatched++;

                    //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
                }
                for(int j = 3; j <= z; j++) {
                    int nlny = (int)(0.5f + ppl.PRECISIONFACTOR*(nlyMass + MassSpecConstants.MASSPROTON)/j);
                    if(nlny < lastTrue && nlny > firstTrue) {
                        numTheroticPeaks++;
                        if(boolMasses[nlny]) numPeaksMatched++;
                    }
                }
            }

            
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }

    public void getNLOnlyChargeNNumPeaksMatched(int z) {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        // y is used for both y and z, ie cterm ions
        // b is used for both b andc, ie nterm ions
        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function 
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        double yneutralloss = 0;
        double bneutralloss = 0;

        for(int i = start; i < end; i++) {

            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            if(bneutralloss == 0 && bd != null) bneutralloss = bd.getNeutralLoss(); 
            if(yneutralloss == 0 && yd != null) yneutralloss = yd.getNeutralLoss(); 

            double bneutralmass = bMass - bneutralloss;
            double yneutralmass = yMass - yneutralloss;

            int int2B = (int)(bneutralmass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yneutralmass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bneutralmass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yneutralmass + MassSpecConstants.MASSPROTON)/2f);
            // doubly charged fragment ions
            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;
 
            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
           // if(boolMasses[int2B++]) numPeaksMatched++;
           // if(boolMasses[int2Y++]) numPeaksMatched++;
            //numTheroticPeaks += 6;
            for(int j = 3; j <= z; j++) {
                int nb = (int)(0.5f + ppl.PRECISIONFACTOR*(bneutralmass + MassSpecConstants.MASSPROTON)/j);
                int ny = (int)(0.5f + ppl.PRECISIONFACTOR*(yneutralmass + MassSpecConstants.MASSPROTON)/j);
                if(nb < lastTrue && nb > firstTrue) { 
                    numTheroticPeaks++;
                    if(boolMasses[nb]) numPeaksMatched++;
                }
                if(ny < lastTrue && ny > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[ny]) numPeaksMatched++;
                }
            }
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }
    public void getChargeNNumPeaksMatched(int z) {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        // y is used for both y and z, ie cterm ions
        // b is used for both b andc, ie nterm ions
        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function 
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;
        for(int i = start; i < end; i++) {


            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);
//System.out.println((char)seq[i] + "\t" + bMass + "\t" + yMass);
            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/2f);
            // doubly charged fragment ions
            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;
 
            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
           // if(boolMasses[int2B++]) numPeaksMatched++;
           // if(boolMasses[int2Y++]) numPeaksMatched++;
            //numTheroticPeaks += 6;
            for(int j = 3; j <= z; j++) {
                int nb = (int)(0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/j);
                int ny = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/j);
                if(nb < lastTrue && nb > firstTrue) { 
                    numTheroticPeaks++;
                    if(boolMasses[nb]) numPeaksMatched++;
                }
                if(ny < lastTrue && ny > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[ny]) numPeaksMatched++;
                }
            }
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }
    public void getNLCharge3NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion
        
        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();

        double yneutralloss = 0;
        double bneutralloss = 0;

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 

            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/2f);

            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;

            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

            if(bneutralloss == 0 && bd != null) bneutralloss = bd.getNeutralLoss(); 
            if(yneutralloss == 0 && yd != null) yneutralloss = yd.getNeutralLoss(); 

            if(bneutralloss != 0) {
                double bneutralmass = bMass - bneutralloss;
                int int2nlB = (int)(bneutralmass*ppl.PRECISIONFACTOR + 0.5f);
                int int3nlB = (int) (0.5f + ppl.PRECISIONFACTOR*(bneutralmass + MassSpecConstants.MASSPROTON)/2f);
                if(int3nlB < lastTrue && int3nlB > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int3nlB]) numPeaksMatched++;
                }
                // singly charged fragment ions
                if(int2nlB < lastTrue && int2nlB > firstTrue) {
                    //numTheroticPeaks += 2;
                    numTheroticPeaks++;
                    if(boolMasses[int2nlB]) numPeaksMatched++;
                }
            }
           if(yneutralloss != 0) { 
                double yneutralmass = yMass - yneutralloss;
                int int2nlY = (int)(yneutralmass*ppl.PRECISIONFACTOR + 0.5f);
                int int3nlY = (int)(0.5f + ppl.PRECISIONFACTOR*(yneutralmass + MassSpecConstants.MASSPROTON)/2f);
                if(int3nlY < lastTrue && int3nlY > firstTrue) {
                    numTheroticPeaks++;
                    if(boolMasses[int3nlY]) numPeaksMatched++;
                }
                if(int2nlY < lastTrue && int2nlY > firstTrue) {
                    //numTheroticPeaks += 2;
                    numTheroticPeaks++;
                    if(boolMasses[int2nlY]) numPeaksMatched++;

                }
           }
        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }
    public void getNLOnlyCharge3NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();

        double yneutralloss = 0;
        double bneutralloss = 0;

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  

            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            if(bneutralloss == 0 && bd != null) bneutralloss = bd.getNeutralLoss(); 
            if(yneutralloss == 0 && yd != null) yneutralloss = yd.getNeutralLoss(); 

            double bneutralmass = bMass - bneutralloss;
            double yneutralmass = yMass - yneutralloss;

            int int2B = (int)(bneutralmass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yneutralmass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bneutralmass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yneutralmass + MassSpecConstants.MASSPROTON)/2f);

            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;

            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }
    public void getCharge3NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        double yMass = getCTermStart(); // overriden function
        double bMass = getNTermStart(); // overriden function
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion
        

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSPROTON)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSPROTON)/2f);

            if(int3B < lastTrue && int3B > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3B]) numPeaksMatched++;

            }
            if(int3Y < lastTrue && int3Y > firstTrue) {
                numTheroticPeaks++;
                if(boolMasses[int3Y]) numPeaksMatched++;
            }
            // singly charged fragment ions
            if(int2B < lastTrue && int2B > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2B]) numPeaksMatched++;
                //if(boolMasses[int2B+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(int2Y < lastTrue && int2Y > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[int2Y]) numPeaksMatched++;

                //if(boolMasses[int2Y+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());

    }

    public void getCharge2NumPeaksMatched(int mam) {
        switch(mam) {
            case 1 : getNLCharge2NumPeaksMatched(); break;// both NL and no NL peaks
            case 2 : getNLOnlyCharge2NumPeaksMatched(); break;
            default : getCharge2NumPeaksMatched(); break;

        }
    }
    public void getCharge3NumPeaksMatched(int mam) {

        switch(mam) {
            case 1 : getNLCharge3NumPeaksMatched(); break;// both NL and no NL peaks
            case 2 : getNLOnlyCharge3NumPeaksMatched(); break;
            default : getCharge3NumPeaksMatched(); break;

        }
    }
    public void getChargeNNumPeaksMatched(int z, int mam) {

        switch(mam) {
            case 1 : getNLChargeNNumPeaksMatched(z); break;// both NL and no NL peaks
            case 2 : getNLOnlyChargeNNumPeaksMatched(z); break;
            default : getChargeNNumPeaksMatched(z); break;

        }
    }
    public void getNLCharge2NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        byte [] seq = getParent().getSequenceAsBytes();
        boolean [] boolMasses = ppl.getBoolMasses();
        MassCalculator mc = ppl.getMassCalculator();
        //fragMasses = mc.getFragMasses(); no big difference
        double yMass = getCTermStart(); // from super class
        double bMass = getNTermStart();  // from super class

        int start = peptide.getStart();
        int end = peptide.getEnd();
        int yIndex = start + end; // index for y ion, i is used for b ion

        double yneutralloss = 0;
        double bneutralloss = 0;

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);
           
            int tempB = (int)(bMass * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)(yMass * ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            if(tempB < lastTrue && tempB > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempB]) numPeaksMatched++;
                //if(boolMasses[tempB+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(tempY < lastTrue && tempY > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempY]) numPeaksMatched++;

                //if(boolMasses[tempY+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

            if(bneutralloss == 0 && bd != null) bneutralloss = bd.getNeutralLoss(); 
            if(yneutralloss == 0 && yd != null) yneutralloss = yd.getNeutralLoss(); 

            if(bneutralloss != 0) {
                int tempnlB = (int)((bMass - bneutralloss) * ppl.PRECISIONFACTOR+ 0.5f);
                if(tempnlB < lastTrue && tempnlB > firstTrue) {
                    //numTheroticPeaks += 2;
                    numTheroticPeaks++;
                    if(boolMasses[tempnlB]) numPeaksMatched++;
                    //if(boolMasses[tempB+ppl.PRECISIONFACTOR]) numPeaksMatched++;
                }
            }
            if(yneutralloss != 0) {
                int tempnlY = (int)((yMass - yneutralloss)*ppl.PRECISIONFACTOR + 0.5f);
                if(tempnlY < lastTrue && tempnlY > firstTrue) {
                    //numTheroticPeaks += 2;
                    numTheroticPeaks++;
                    if(boolMasses[tempnlY]) numPeaksMatched++;

                    //if(boolMasses[tempY+ppl.PRECISIONFACTOR]) numPeaksMatched++;
                }
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());
    }
    public void getNLOnlyCharge2NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        byte [] seq = getParent().getSequenceAsBytes();
        boolean [] boolMasses = ppl.getBoolMasses();
        MassCalculator mc = ppl.getMassCalculator();
        //fragMasses = mc.getFragMasses(); no big difference
        double yMass = getCTermStart(); // from super class
        double bMass = getNTermStart();  // from super class

        int start = peptide.getStart();
        int end = peptide.getEnd();
        int yIndex = start + end; // index for y ion, i is used for b ion

        double bneutralloss = 0; // neutral loss b
        double yneutralloss = 0; // neutral loss y

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 

            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);

            yMass += mc.getFragmentMass(seq[yIndex-i]);

            if(bneutralloss == 0 && bd != null) bneutralloss = bd.getNeutralLoss(); 
            if(yneutralloss == 0 && yd != null) yneutralloss = yd.getNeutralLoss(); 

            int tempB = (int)((bMass - bneutralloss) * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)((yMass - yneutralloss) * ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            if(tempB < lastTrue && tempB > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempB]) numPeaksMatched++;
                //if(boolMasses[tempB+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(tempY < lastTrue && tempY > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempY]) numPeaksMatched++;

                //if(boolMasses[tempY+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());
    }
    public void getCharge2NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        int firstTrue = ppl.getFirstTrue() - 1;
        int lastTrue = ppl.getLastTrue() + 1;

        byte [] seq = getParent().getSequenceAsBytes();
        boolean [] boolMasses = ppl.getBoolMasses();
        MassCalculator mc = ppl.getMassCalculator();
        //fragMasses = mc.getFragMasses(); no big difference
        double yMass = getCTermStart(); // from super class
        double bMass = getNTermStart();  // from super class

        int start = peptide.getStart();
        int end = peptide.getEnd();
        int yIndex = start + end; // index for y ion, i is used for b ion

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);

            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int tempB = (int)(bMass * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            if(tempB < lastTrue && tempB > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempB]) numPeaksMatched++;
                //if(boolMasses[tempB+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }
            if(tempY < lastTrue && tempY > firstTrue) {
                //numTheroticPeaks += 2;
                numTheroticPeaks++;
                if(boolMasses[tempY]) numPeaksMatched++;

                //if(boolMasses[tempY+ppl.PRECISIONFACTOR]) numPeaksMatched++;
            }

        }
        setNumPeaksMatched(numPeaksMatched, numTheroticPeaks, ppl.getPTrue());
    }

    public int [] getChargeNTheorMasses(int z, int mam) { 

        switch(mam) {
            case 1 : return getNLChargeNTheorMasses(z); // both NL and no NL peaks
            case 2 : return getNLOnlyChargeNTheorMasses(z);
            default : return getChargeNTheorMasses(z);

        }

    }
    public int[] getCharge3TheorMasses(int mam) { 

        switch(mam) {
            case 1 : return getNLCharge3TheorMasses(); // both NL and no NL peaks
            case 2 : return getNLOnlyCharge3TheorMasses();
            default : return getCharge3TheorMasses();

        }

    }
    public int [] getCharge2TheorMasses(int mam) {
        switch(mam) {
            case 1 : return getNLCharge2TheorMasses(); // both NL and no NL peaks
            case 2 : return getNLOnlyCharge2TheorMasses();
            default : return getCharge2TheorMasses();

        }
    }
    public int [] getNLCharge2TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion

        double bneutralloss = 0; // neutral loss b
        double yneutralloss = 0; // neutral loss y

        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {
            int yi = yIndex -i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }
            // specific for diff mod
            DiffMod bm = diffMods[i];
            DiffMod ym = diffMods[yi];

            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();


            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            if(bneutralloss == 0 && bm != null) bneutralloss = bm.getDbinwidthNeutralLoss(); 
            if(yneutralloss == 0 && ym != null) yneutralloss = ym.getDbinwidthNeutralLoss(); 

            if(bneutralloss != 0) processSinglyChargedBIon(theorMass, bMass - bneutralloss);
            if(yneutralloss != 0) processSinglyChargedYIon(theorMass, yMass - yneutralloss);
        }

        return theorMass;
    }
    public int [] getNLOnlyCharge2TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
       
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion

        double bneutralloss = 0; // neutral loss b
        double yneutralloss = 0; // neutral loss y
        
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {
            int yi = yIndex -i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }
            // specific for diff mod
            DiffMod bm = diffMods[i];
            DiffMod ym = diffMods[yi];

            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();
            //yMass -= ym == null? 0 : ym.getDbinwidthNeutralLoss();

            if(bneutralloss == 0 && bm != null) bneutralloss = bm.getDbinwidthNeutralLoss(); 
            if(yneutralloss == 0 && ym != null) yneutralloss = ym.getDbinwidthNeutralLoss(); 

            processSinglyChargedBIon(theorMass, bMass - bneutralloss);
            processSinglyChargedYIon(theorMass, yMass - yneutralloss);
        }

        return theorMass;
    }
    public int [] getCharge2TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
       
        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {
            int yi = yIndex -i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }
            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);
        }

        return theorMass;
    }
    // both neutral loss and no neutral loss fragement peaks
    public int[] getNLCharge3TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion

        double yneutralloss = 0;
        double bneutralloss = 0;

        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex -i;

            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];
            
            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod
            DiffMod bm = diffMods[i];
            DiffMod ym = diffMods[yi];

            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);
           
            if(bneutralloss == 0 && bm != null) bneutralloss = bm.getDbinwidthNeutralLoss(); 
            if(yneutralloss == 0 && ym != null) yneutralloss = ym.getDbinwidthNeutralLoss(); 

            if(bneutralloss != 0) {             
                processSinglyChargedBIon(theorMass, bMass - bneutralloss);
                processDoublyChargedBIon(theorMass, bMass - bneutralloss);
            }
            if(yneutralloss != 0) {
                processSinglyChargedYIon(theorMass, yMass - yneutralloss);
                processDoublyChargedYIon(theorMass, yMass - yneutralloss);
            }
        }
        return theorMass;
    }
    public int[] getNLOnlyCharge3TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();

        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion

        double bneutralloss = 0; // neutral loss b
        double yneutralloss = 0; // neutral loss y

        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex -i;

            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];
     
            
            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod
            DiffMod bm = diffMods[i];
            DiffMod ym = diffMods[yi];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            if(bneutralloss == 0 && bm != null) bneutralloss = bm.getDbinwidthNeutralLoss(); 
            if(yneutralloss == 0 && ym != null) yneutralloss = ym.getDbinwidthNeutralLoss(); 

            processSinglyChargedBIon(theorMass, bMass - bneutralloss);
            processSinglyChargedYIon(theorMass, yMass - yneutralloss);

            processDoublyChargedBIon(theorMass, bMass - bneutralloss);
            processDoublyChargedYIon(theorMass, yMass - yneutralloss);
        }
        return theorMass;
    }
    public int[] getCharge3TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();
        
        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex -i;

            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];
     
            
            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);
        }
        return theorMass;
    }

    public int [] getNLOnlyChargeNTheorMasses(int z) {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;


        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion

        double yneutralloss = 0;
        double bneutralloss = 0;

        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex-i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod

            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            if(bneutralloss == 0 && bm != null) bneutralloss = bm.getDbinwidthNeutralLoss(); 
            if(yneutralloss == 0 && ym != null) yneutralloss = ym.getDbinwidthNeutralLoss(); 

            processSinglyChargedBIon(theorMass, bMass - bneutralloss);
            processSinglyChargedYIon(theorMass, yMass - yneutralloss);

            processDoublyChargedBIon(theorMass, bMass - bneutralloss);
            processDoublyChargedYIon(theorMass, yMass - yneutralloss);

            // for highly charged fragment ions
            processHighlyChargedIons(theorMass, bMass - bneutralloss, yMass - yneutralloss, z);
        }
        return theorMass;
    }
    public int [] getNLChargeNTheorMasses(int z) {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;

        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion

        double yneutralloss = 0;
        double bneutralloss = 0;

        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex-i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod
            DiffMod bm = diffMods[i];
            DiffMod ym = diffMods[yi];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();


            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);

            // for highly charged fragment ions
            processHighlyChargedIons(theorMass, bMass, yMass, z);

            if(bneutralloss == 0 && bm != null) bneutralloss = bm.getDbinwidthNeutralLoss(); 
            if(yneutralloss == 0 && ym != null) yneutralloss = ym.getDbinwidthNeutralLoss(); 
            double nlbMass = bMass - bneutralloss;
            double nlyMass = yMass - yneutralloss;

            if(bneutralloss != 0) {
                processSinglyChargedBIon(theorMass, nlbMass);
                processDoublyChargedBIon(theorMass, nlbMass);
            }
            if(yneutralloss != 0) {
                processSinglyChargedYIon(theorMass, nlyMass);
                processDoublyChargedYIon(theorMass, nlyMass);
            }

            // for highly charged fragment ions
            if(bneutralloss != 0 || yneutralloss != 0) {
                processHighlyChargedIons(theorMass, nlbMass, nlyMass, z);
            }
        }
        return theorMass;
    }
    public int [] getChargeNTheorMasses(int z) {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        byte [] seq = getSequence().getBytes();
        double bMass = getNTermDbinwidthStart() + 0.5f;
        double yMass = getCTermDbinwidthStart() + 0.5f;


        int leftshift1 = ppl.getSearchParams().getLeftShift1();
        int leftshift2 = ppl.getSearchParams().getLeftShift2();
        int leftshift3 = ppl.getSearchParams().getLeftShift3();

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {

            int yi = yIndex-i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            if(i == leftshift1 || i == leftshift2 || i == leftshift3) {
                bMass -= 1;
                yMass -= 1;
            }

            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();


            processSinglyChargedBIon(theorMass, bMass);
            processSinglyChargedYIon(theorMass, yMass);

            processDoublyChargedBIon(theorMass, bMass);
            processDoublyChargedYIon(theorMass, yMass);

            // for highly charged fragment ions
            processHighlyChargedIons(theorMass, bMass, yMass, z);
        }
        return theorMass;
    }
    public boolean equals(ModifiedPeptideHit m) {
        return super.equals(m) && mods.equals(m.mods);
    }


}
