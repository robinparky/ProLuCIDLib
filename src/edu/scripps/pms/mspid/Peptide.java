/**
 * @file Peptide.java
 * This is the source file for edu.scripps.pms.util.seq.Peptide
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:41:49 $
 */



package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;

import java.util.List;

public class Peptide {

    // description line of this Peptide sequence
    // the sequence string of this Peptide
    
    protected int start;
    protected int end;
    protected Fasta parent;
    //private int chargeState;

    public Peptide(Fasta parent, int start, int end) {
        this.parent = parent;
        this.start = start;
        this.end = end;
    }  

    public boolean equals(Peptide p) {

        if(p != null) {
            return parent == p.parent && start == p.start && end == p.end; 
        }
        return false;
    }
    public Fasta getParent() {
        return parent;
    }
    public String getExtraExtendedSequence() {
        // peptide sequence with 3 more residues before and after the '.'
        // 8 more bytes for the heading and tailing aa and '.'
        int len = getLength() + 8;
        byte [] seq = new byte[len];
        seq[0] = '-';
        seq[1] = '-';
        seq[2] = '-';
        seq[len-1] = '-';
        seq[len-2] = '-';
        seq[len-3] = '-';

        switch(start){
            case 0: 
                break;
            case 1: 
                seq[2] = parent.byteAt(start-1);
                break;
            case 2: 
                seq[1] = parent.byteAt(start-2);
                seq[2] = parent.byteAt(start-1);
                break;
            default: 
                seq[0] = parent.byteAt(start-3);
                seq[1] = parent.byteAt(start-2);
                seq[2] = parent.byteAt(start-1);
                break;

        }
        seq[3] = '.';

        int lastIndex = parent.getLength() - 1;
        int diff = lastIndex - end;


        switch(diff){
            case 0: 
                break;
            case 1: 
                seq[len-3] = parent.byteAt(end+1);
                break;
            case 2: 
                seq[len-3] = parent.byteAt(end+1);
                seq[len-2] = parent.byteAt(end+2);
                break;
            default: 
                seq[len-3] = parent.byteAt(end+1);
                seq[len-2] = parent.byteAt(end+2);
                seq[len-1] = parent.byteAt(end+3);
                break;

        }
        seq[len-4] = '.';
        int index = 4;
        for (int i = start; i < end+1; i++) {
            seq[index++] = parent.byteAt(i);
        }

        return new String(seq);
    }
    public String getExtendedSequence() {

        // 4 more bytes for the heading and tailing aa and '.'
        int len = getLength() + 4;
        byte [] seq = new byte[len];
        if (start == 0) {
            seq[0] = '-';
        } else {
            seq[0] = parent.byteAt(start-1);
        }
        seq[1] = '.';

        int lastIndex = parent.getLength() - 1;
        if (end == lastIndex) {
            seq[len-1] = '-';
        } else {
            seq[len-1] = parent.byteAt(end+1);
        }
        seq[len-2] = '.';
        int index = 2;
        for (int i = start; i < end+1; i++) {
            seq[index++] = parent.byteAt(i);
        }

        return new String(seq);
    }
    public String getSequence() {

        try {
            new String(parent.getSequenceAsBytes(), start, getLength());
        } catch(Exception e) {
            e.printStackTrace();
        }
        return new String(parent.getSequenceAsBytes(), start, getLength());
    }
    public int getStart() {
        return start;
    }
    public String getDefline() {
        return parent.getDefline();
    }
    public int getEnd() {
        return end;
    }

    public int getLength() {
        return end - start + 1;
    }
    public String getSequestLikeAccession() {
        return parent.getSequestLikeAccession();
    }
    public List<String> getDefList() {
        return parent.getDefList();
    }
    public String getAccession() {
        return parent.getAccession();
    }
}
