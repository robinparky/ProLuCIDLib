package edu.scripps.pms.util.enzyme;

import edu.scripps.pms.mspid.PeptideHit;
import edu.scripps.pms.util.seq.Fasta;

/**
 * @author  Tao Xu
 * @version $Id
 */

public class Protease {
    private static int NUMCHARS = 256;

    private String name;
    // true for cut at C terminus and false for N terminus
    private boolean isC;
    // list of cleavage sites
    private boolean [] residues = new boolean[NUMCHARS];

    
    public Protease(String name) {
        this.name = name;
    }
    // 
    public void setType(boolean isC) {
        this.isC = isC;
    }
    public boolean getType() {
        return isC;
    }
    public void addCleavageSite(char residue) {
        residues[residue] = true;
    }
    public int checkEnzymeSpecificity(PeptideHit p) {
        byte [] seq = p.getParent().getSequenceAsBytes();
        return checkLeft(seq, p.getStart()) + checkRight(seq, p.getEnd());
    }
    public int checkEnzymeSpecificity(Fasta sequence, int start, int end) {
        byte [] seq = sequence.getSequenceAsBytes();
        return checkLeft(seq, start) + checkRight(seq, end);
    }
    public int checkEnzymeSpecificityStrict(PeptideHit p) {
        byte [] seq = p.getParent().getSequenceAsBytes();
        return checkLeftStrict(seq, p.getStart()) + checkRight(seq, p.getEnd());
    }
    public int checkEnzymeSpecificityStrict(Fasta sequence, int start, int end) {
        byte [] seq = sequence.getSequenceAsBytes();
        return checkLeftStrict(seq, start) + checkRight(seq, end);
    }
    public int checkEnzymeSpecificityStrict(byte [] seq, int start, int end) {
        
        return checkLeftStrict(seq, start) + checkRight(seq, end);
    }
    public int getNumInternalMissCleavage(Fasta sequence, int start, int end) {
        byte [] seq = sequence.getSequenceAsBytes();
        return getNumInternalMissCleavage(seq, start, end);
    }
    public int getNumInternalMissCleavage(byte [] seq, int start, int end) {
        int result = 0;
        if(isC) {
            end--;
        } else {
            start++;
        }       
  
        for(int i = start; i <= end; i++) {
            if(residues[seq[i]]) {
                result++;
            }
        }
        return result;     
    }
    public boolean isDigestable(char c) {
        return residues[c];
    }
    public String getResidues() {
        StringBuffer sb = new StringBuffer(10);
        for(int i = 0; i < residues.length; i++) {
            if(residues[i]) {
                sb.append((char)i);
            }
        }
        return sb.toString();
    }
    public int checkLeftStrict(byte [] seq, int index) {
        if(index == 0) { 
            return 1;
        }
        if(isC) {
            return residues[seq[--index]]? 1 : 0;
        } else {
            return residues[seq[index]]? 1 : 0;
        }

    }
    public int checkLeft(byte [] seq, int index) {
        //if(index <= 30) { // for N-term truncation, like signal peptide etc
        if(index <= 1) { // for N-term truncation, like signal peptide etc
            return 1;
        }
        if(String.valueOf((char)(seq[index-1])).equals("-")){
            return 1;
        }
        if(isC) {
            return residues[seq[--index]]? 1 : 0;
        } else {

            return residues[seq[index]]? 1 : 0;
        }

    }
    public int checkRight(byte [] seq, int index) {
        // should we also consder C-term truncation?
        if(String.valueOf((char)(seq[index+1])).equals("-")){
            return 1;
        }
        
        
        
        if(index >= seq.length - 1) { 
            return 1; // last residue
        }
      
        // now not last residue
        if(isC) {
            return residues[seq[index]]? 1 : 0;
        } else {
            return residues[seq[++index]]? 1 : 0;
        }
    }
    public String getName() {
        return name;
    }

    public static void main(String[] args) {
       Protease p = new Protease("Trypsin");
       p.setType(true);
              
       p.addCleavageSite('R');
       p.addCleavageSite('K');
        
       
    }
}
