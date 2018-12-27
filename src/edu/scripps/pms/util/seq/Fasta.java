/**
 * @file Fasta.java
 * This is the source file for edu.scripps.pms.util.spectrum.Fasta
 * @author Tao Xu
 * @author Robin Park
 * @date $Date: 2014/08/15 18:05:41 $
 */



package edu.scripps.pms.util.seq;

import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

public class Fasta implements Comparable<Fasta> {

    // description line of this Fasta sequence
    //protected String defline;
    protected byte [] def;
    // the sequence string of this Fasta
    protected byte [] sequence;
    //private String seq;
    //private String accession = null;
    //private String sequestLikeAccession = null;
    private double mPlusH = 0;
    private static final Pattern pattern = Pattern.compile("(.*)\\d");
    protected List<String> defList = new ArrayList<String>();

    public Fasta(String sequence) {
        sequence = sequence.toUpperCase();
        try {
            this.sequence = sequence.getBytes("US-ASCII");
        } catch (UnsupportedEncodingException e) {
            System.err.println("Unknow charset");
            System.exit(1);
        }
    }

    public void addDefList(String def) {
        this.defList.add(def);
    }

    public Fasta(String defline, String sequence) {
        //this.defline = defline;
//System.out.println("defline: " + this.defline);
        sequence = sequence.toUpperCase();
    //    seq = sequence;
        try {
         
            def = defline.getBytes("US-ASCII");
            this.sequence = sequence.getBytes("US-ASCII");
        } catch (UnsupportedEncodingException e) {
            System.err.println("Unknow charset");
            System.exit(1);
        }
//System.out.println("hahaha, defline: " + this.defline);
    }
    public void setMPlusH(double mh) {
        mPlusH = mh;
    }
    public double getMPlusH() {
        return mPlusH;
    }
    //public double setMPlusH() {
    //    return mPlusH;
    //}
    // to order the sequence from long to short
    public int compareTo(Fasta f) {
        if(f == null) return -1;
        return f.getLength() - getLength(); 
    }

    // return null if Gene_Symbol is not included in defline
    public String getGeneName() {
        String defline = getDefline();
        String genenamestring = defline.indexOf("Gene_Symbol=") > -1? "Gene_Symbol=" : "GN=";
        if(defline.indexOf(genenamestring) == -1) return null; // neither Gene_Symbol= nor GN= is found.
        //String [] arr = defline.split("Gene_Symbol=");    
        String[] arr = defline.split(genenamestring);
        if(arr.length > 1) {
            String geneele = arr[1];
            String genename = geneele.split("\\s+")[0];
//System.out.println(geneele + "\t" + genename + "\t" + defline);
 //System.out.println("gene element: " + geneele + "\t" + genename + "Haha");
            return genename;

        } else {

//System.out.println("no gene nmae: "  +  defline);
            return null;
        }
    }

    public String getAccessionWithNoVersion() {
        String accession = getAccession();
        int index = accession.indexOf(".");
        if(index == -1) {
            return accession;
        } else {
            return getAccession().substring(0, index);
        }
    }
    public String getSequence() {
        return new String(sequence);
        //return new String(sequence);
    }
    public byte [] getSequenceAsBytes() {
        return sequence;
    }
    public String getOriginalDefline() {
        
        //return defline.substring(1,defline.length());
        return new String(def);
    }
    public String getDefline() {
        
        return getOriginalDefline().substring(1,def.length);
    }

    public String getDescription() {
        String origdefline = getOriginalDefline();
        int space = origdefline.indexOf(" ");
        int tab = origdefline.indexOf("\t");
        
        int index = 1;
        if(tab > 0) {
            if(space > 0) {
                index = tab > space? space : tab;
            } else {
                index = tab + 1;
            }
        } else {
            if(space > 0) {
                index = space + 1;
            }
        }
        
        return origdefline.substring(index, def.length);
    }
    public byte byteAt(int index) {
        return sequence[index];
    }
    public int getLength() {
        return sequence.length;
    }


    // get accession without version
    public String getAccession() {
//System.out.println("defline: " + defline);
        return getAccession(getOriginalDefline().substring(1));

    }

    // get the first refseq accession number in the defline
    public String getEnsemblsAccession() {
        String defline = getOriginalDefline();
        String locus = defline.split(" ")[0];
        String[] arr = locus.split("\\|");
        for(int i = 0; i < arr.length; i++) {
            String s = arr[i];
            if(s.startsWith("ENSEMBL")) {
                String refacs = s.split(":")[1];
                
                return (refacs.split(";")[0]);
            }
        }
        return null;
    }
    // get the first refseq accession number in the defline
    public String getTremblAccession() {
        String defline = getOriginalDefline();
        String locus = defline.split(" ")[0];
        String[] arr = locus.split("\\|");
        for(int i = 0; i < arr.length; i++) {
            String s = arr[i];
            if(s.startsWith("TREMBL")) {
                String refacs = s.split(":")[1];
                
                return (refacs.split(";")[0]);
            }
        }

        return null;
    }
    // get the first refseq accession number in the defline
    public String getVegaAccession() {
        String defline = getOriginalDefline();
        String locus = defline.split(" ")[0];
        String[] arr = locus.split("\\|");
        for(int i = 0; i < arr.length; i++) {
            String s = arr[i];
            if(s.startsWith("VEGA")) {
                String refacs = s.split(":")[1];
                
                return (refacs.split(";")[0]);
            }
        }

        return null;
    }
    // get the first refseq accession number in the defline
    public String getSwissprotAccession() {
        String defline = getOriginalDefline();
        String locus = defline.split(" ")[0];
        String[] arr = locus.split("\\|");
        for(int i = 0; i < arr.length; i++) {
            String s = arr[i];
            if(s.startsWith("SWISS-PROT")) {
                String refacs = s.split(":")[1];
                
                return (refacs.split(";")[0]);
            }
        }

        return null;
    }
    // get the first refseq accession number in the defline
    public String getRefSeqAccession() {
        StringBuffer sb = new StringBuffer();
        String defline = getOriginalDefline();
        String locus = defline.split(" ")[0];
        String[] arr = locus.split("\\|");
        for(int i = 0; i < arr.length; i++) {
            String s = arr[i];
            if(s.startsWith("REFSEQ")) {
                String refacs = s.split(":")[1];
                
                return (refacs.split(";")[0]);
            }
        }

        return null;
    }
    // get the refseq accession numbers (separated with a space char if multiple) from the IPI databases
    public String getAllRefSeqAccessions() {
        StringBuffer sb = new StringBuffer();
        String defline = getOriginalDefline();
        String locus = defline.split(" ")[0];
        String[] arr = locus.split("\\|");
        for(int i = 0; i < arr.length; i++) {
            String s = arr[i];
            if(s.startsWith("REFSEQ")) {
                String refacs = s.split(":")[1].replace(';', ' ');
                sb.append(refacs);
            }
        }

        return sb.toString();
    }
   
    public String getSequestLikeAccession() {
        return getSequestLikeAccession(getOriginalDefline().substring(1));
    }
    public boolean isReversed() {
        return getAccession().startsWith("Re");
    }
    public static String getSequestLikeAccession(String acc) {
        String[] arr = acc.split("\t");
        String[] arr1 = arr[0].split(" ");
        String newacc = arr1[0];
        if(newacc != null && newacc.length() > 40) {
            newacc = newacc.substring(0, 40);
        }
        return newacc;        

    }

    public String getLongAccession() {
        return getLongAccession(getOriginalDefline().substring(1));
    }
    public static String getLongAccession(String acc) {
        String[] arr = acc.split("\t");
        String[] arr1 = arr[0].split(" ");
        String newacc = arr1[0];
        return newacc;        
    }

    public static String getAccession(String accession)
    {
        //NCBI, IPI, or others such as UNIT_PROT, SGD, NCI
//        accession = getDefline().substring( getDefline().indexOf('>')+1 );
        //accession = getDefline();

        //There are many corruptted sqt file.  Ignore it.
        try
        {
            if( accession.startsWith("gi") && accession.contains("|") ) //NCBI
            {
                String[] arr = accession.split("\\|+");

		if( arr.length>=4 && ("gb".equals(arr[2]) || "ref".equals(arr[2]) || "emb".equals(arr[2]) || "dbj".equals(arr[2]) || "prf".equals(arr[2]) ||"sp".equals(arr[2])) || "tpd".equals(arr[2]) ||"tpg".equals(arr[2]) ||"tpe".equals(arr[2]) )
		    accession = arr[3];
                else
                {
                    arr = accession.split(" ");
                    accession = arr[0];
                }

                //Accession # should end with digit.  If accession # does not end with digit,
                //grap next string (We assume this next one ends with digit.)
		/*
                if( pattern.matcher(arr[3]).matches() )
                    accession = arr[3];
                else
                    accession = arr[4].substring(0, arr[4].indexOf(" "));
		*/

            }
            else if( accession.startsWith("sp") || accession.startsWith("tr")) //IPI
            {
                String[] arr = accession.split("\\|");
                accession = arr[1];
            }
            else if( accession.startsWith("IPI") ) //IPI
            {
                String arr[] = accession.split("\\|");
                String subArr[] = arr[0].split(":");

                if(subArr.length>1)
                    accession = subArr[1];
                else
                    accession = subArr[0];
            }
            else if( accession.startsWith("Re") || accession.startsWith("contam") || accession.startsWith("Contam")) //Reverse database
            {
                int space = accession.indexOf(" ");
                int tab = accession.indexOf("\t");

                if(space<0) space = 40;
                if(tab<0) tab = 40;

                int index = (tab>space)?space:tab;

                int end;

                if(index<=0 || index>=40) //no space
                {
                    int length = accession.length();
                    end = (length>40)?40:length;
                }
                else  //cut by the first space
                    end = index;

                accession = accession.substring(0, end);
            }
            else //UNIT_PROT, NCI or SGD
            {
                int spaceIndex = accession.indexOf(" ");
                int tabIndex = accession.indexOf("\t");

                if(spaceIndex>0)
                {

                    if(tabIndex>0 && spaceIndex>tabIndex)
                        accession = accession.substring(0, tabIndex);
                    else
                        accession = accession.substring(0, spaceIndex);
                } else { // no space
                    
                    if(tabIndex>0)
                        accession = accession.substring(0, tabIndex);
                }
            }
        }
        catch(Exception e)
        {
            //System.out.println("No Correct Accession found, but this will be handled by MSP system." + accession + " " +  e);

            int i = accession.indexOf(" ");
            if(i<0)
                return accession;
            else
                return accession.substring(0, i);

        }

        return accession;
    }

    public static void main(String args[]) throws Exception
    {
        String seq = getSequestLikeAccession("IPI:IPI00004273.6|SWISS-PROT:P49756-1|TR");
        System.out.println(seq);

    }

    public List<String> getDefList() {
        return defList;
    }

    public void setDefList(List<String> defList) {
        this.defList = defList;
    }
}
