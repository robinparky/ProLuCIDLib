package edu.scripps.dia;

import edu.scripps.dia.util.IndexUtil;
import edu.scripps.dia.util.IndexedFile;
import gnu.trove.TFloatArrayList;
import gnu.trove.TLongHashSet;
import org.jdom.JDOMException;
import org.sqlite.SQLiteConfig;

import java.io.*;
import java.nio.ByteBuffer;
import java.sql.*;
import java.util.*;

import static edu.scripps.pms.util.seq.Fasta.getSequestLikeAccession;

/**
 * Created by Titus Jung titusj@scripps.edu on 3/26/18.
 */
public class LibraryIndexer {

    private PreparedStatement addEntryStatement;
    private PreparedStatement addSpectraStatement;
    private PreparedStatement getEntryByMassRangeStatement;
    private PreparedStatement getEntrySequenceCSStatement;
    private PreparedStatement updateEntry;
    private PreparedStatement getSpectraStatement;
    private PreparedStatement incrementCopies;
    private PreparedStatement getSequenceCSIDScoreStatement;
    private PreparedStatement getScoreStatement;
    protected PreparedStatement updateSpectraStatement;
    protected Connection con;
    private long peptideID =0;
    private long proteinID = 0;
    private Map<String, IndexedFile> indexedMap;
    private PreparedStatement getProteinIDStatement;
    private PreparedStatement getPeptideCountStatement;
    private PreparedStatement getProteinCountStatement;
    private PreparedStatement addProteinEntryStatement;
    private PreparedStatement addProteinToPeptideStatement;
    private PreparedStatement getSpectraMassStatement;
    private PreparedStatement getProteinIDFromLocusStatement;
    private PreparedStatement getPeptideProteinTable;
    private PreparedStatement retrievePeptideIDs;
    private PreparedStatement retreiveSequenceKeys;
    private String dbAddress="";
    private int descriptiveIndex=-1;
    private int locusIndex=-1;
    private int seqIndex=-1;
    private int fileNameIndex=-1;
    private int theoreticalMassIndex=-1;
    private int xcorrIndex = -1;
    private int deltaCNIndex = -1;
    private Map<String, Long>           seqKeyIDMap = new HashMap<>();
    private Map<String, Float>          seqKeyScoreMap = new HashMap<>();
    private Map<String, TLongHashSet>   seqKeyProteinSetMap = new HashMap<>();
    private Map<String, Long>           locusProteinIDMap = new HashMap<>();
    public static final float NO_FLOAT_VALUE = -2000;
    public static final String CSV_HEADER =
            "H\tPeptideID\tFilePath\tSequence\tChargeState\tRetentionTime\tMass\tNumberOfPeaks\tMassList\tIntensityList";


    public final static int MZ_KEY_SHIFT = 1000;

    public final static String PROTEIN_TABLE_NAME = "ProteinTable";
    public final static String PEPTIDE_PROTEIN_INDEX_TABLE = "PeptideProteinIndexTable";
    public static final String proteinSchema =
            "(proteinID," +
                    "Accession," +
                    "Description)";
    public static final String proteinCreationSchema =
            "(proteinID INTEGER," +
                    "Accession TEXT," +
                    "Description TEXT)";

    public static final String peptideProteinIndexTableSchema =
            "(peptideID, " +
                    "proteinID)";

    public static final String peptideProteinIndexTableCreationSchema =
            "(peptideID  INTEGER, " +
                    "proteinID INTEGER)";

    public static final String peptideSchema
            = " (peptideID," + // 1
            "peptideSeq," + // 2
            "massKey," +// 3
            "precursorMZ," +// 4
            "chargeState," +// 5
            "copies," +// 6
            "numPeaks," +// 7
            "retentionTime," +// 8
            "startTime," +// 9
            "endTime," +// 10
            "fileName," +// 11
            "searchScore," +// 12
            "searchScoreType," +// 13
            "deltaCN," +// 14
            "sequenceCS," +
            "scan) ";// 15




    public static final String peptideSchemaCreation
            = " (peptideID INTEGER," +
            "peptideSeq TEXT," +
            "massKey INTEGER," +
            "precursorMZ REAL," +
            "chargeState INTEGER," +
            "copies INTEGER," +
            "numPeaks INTEGER," +
            "retentionTime REAL," +
            "startTime REAL," +
            "endTime REAL," +
            "fileName TEXT," +
            "searchScore REAL," +
            "searchScoreType INTEGER," +
            "deltaCN REAL," +
            "sequenceCS TEXT," +
            "scan INTEGER) ";

    public static final String PEPTIDE_TABLE_NAME = "PeptideTable";
    public static final String SPECTRA_TABLE_NAME = "SpectraTable";
    public static final String spectraShema = " (peptideID,peakMZ,peakIntensity,massKey)";
    public static final String spectraShemaCreation
            = " (peptideID INTEGER," +
            "peakMZ BLOB," +
            "peakIntensity BLOB," +
            "massKey Integer)";






    public LibraryIndexer(String address) throws SQLException {
        dbAddress = address;
        init();
    }

    public LibraryIndexer() {
        seqKeyScoreMap = new HashMap<>();
        seqKeyIDMap = new HashMap<>();
        peptideID = 0;
    }

    protected void init() throws SQLException {
        SQLiteConfig config = new SQLiteConfig();
        final int indexFactor = 1;
        final int cacheSize = 100000 / indexFactor;
        final int pageSize = 4096;
        String tempDir = System.getProperty("java.io.tmpdir");
        config.setTempStoreDirectory(tempDir);
        //optimize for multiple connections that can share data structures
        config.setSharedCache(true);
        config.setCacheSize(cacheSize);
        config.setPageSize(pageSize);
        config.setJournalMode(SQLiteConfig.JournalMode.OFF);
        config.enableEmptyResultCallBacks(false);
        config.enableCountChanges(false);
        config.enableFullSync(false);
        config.enableRecursiveTriggers(false);
        config.setLockingMode(SQLiteConfig.LockingMode.NORMAL);
        config.setSynchronous(SQLiteConfig.SynchronousMode.OFF); //TODO may be dangerous on some systems to have off
        con = DriverManager.getConnection("jdbc:sqlite:" + dbAddress, config.toProperties());

        config.setMaxPageCount(10 * 1000 * 1000 * 1000 / pageSize);

        con.prepareStatement("CREATE TABLE IF NOT EXISTS "+ PEPTIDE_TABLE_NAME +" "+
                peptideSchemaCreation).execute();
        con.prepareStatement("CREATE TABLE IF NOT EXISTS "+PROTEIN_TABLE_NAME+" "+
                proteinCreationSchema).execute();
        con.prepareStatement("CREATE TABLE IF NOT EXISTS "+ SPECTRA_TABLE_NAME
                + " "+spectraShemaCreation).execute();
        con.prepareStatement("CREATE TABLE IF NOT EXISTS "+PEPTIDE_PROTEIN_INDEX_TABLE
                + " "+peptideProteinIndexTableCreationSchema).execute();
        con.prepareStatement("CREATE INDEX IF NOT EXISTS massKey_index_dsc ON "
                + PEPTIDE_TABLE_NAME +" (massKey DESC);").execute();
        con.prepareStatement("CREATE UNIQUE INDEX IF NOT EXISTS sequenceCSKey_index_dsc ON "
                + PEPTIDE_TABLE_NAME +" (sequenceCS DESC);").execute();
        con.prepareStatement("CREATE INDEX IF NOT EXISTS id_index_dsc ON "
                + SPECTRA_TABLE_NAME +" (peptideID DESC);").execute();
        con.prepareStatement("CREATE INDEX IF NOT EXISTS massKey_index_dsc ON "
                + SPECTRA_TABLE_NAME +" (massKey DESC);").execute();
        con.prepareStatement("CREATE INDEX IF NOT EXISTS peptideID_index_dsc ON "
                +  PEPTIDE_PROTEIN_INDEX_TABLE+" (peptideID DESC);").execute();
        con.prepareStatement("CREATE UNIQUE INDEX IF NOT EXISTS proteinID_index_dsc ON "
                +  PROTEIN_TABLE_NAME+" (proteinID DESC);").execute();
        con.prepareStatement("CREATE UNIQUE INDEX IF NOT EXISTS Accession_index_dsc ON "
                +  PROTEIN_TABLE_NAME+" (Accession DESC);").execute();

        addEntryStatement
                = con.prepareStatement("INSERT INTO " + getPeptideTableName() + peptideSchema
                + "VALUES (?, ?, ?, ?, ?, ?, ?, ?,?, ?,?, ?, ?,?,?,?);") ;
        addSpectraStatement
                = con.prepareStatement("INSERT INTO " + SPECTRA_TABLE_NAME + spectraShema
                + "VALUES (?, ?, ?,?);") ;
        getEntryByMassRangeStatement = con.prepareStatement(
                "SELECT * "
                        + "FROM " + PEPTIDE_TABLE_NAME + " "
                        + "INNER JOIN "+ PEPTIDE_PROTEIN_INDEX_TABLE
                        + " ON " +PEPTIDE_TABLE_NAME+".peptideID="+PEPTIDE_PROTEIN_INDEX_TABLE+".peptideID "
                        + "INNER JOIN "+PROTEIN_TABLE_NAME
                        +" ON "+PEPTIDE_PROTEIN_INDEX_TABLE+".proteinID ="+PROTEIN_TABLE_NAME+".proteinID "
                        + "WHERE massKey BETWEEN ? AND ?" +
                        " ORDER BY peptideID;")
        ;
        getEntrySequenceCSStatement= con.prepareStatement(
                "SELECT * "
                        + "FROM " + PEPTIDE_TABLE_NAME + " "
                        + "WHERE sequenceCS = ?;");
        updateEntry = con.prepareStatement(
                "UPDATE "+ PEPTIDE_TABLE_NAME +
                        " SET massKey =?," +
                        "precursorMZ =?," +
                        "numPeaks =?," +
                 //       "copies =?,"+
                        "retentionTime =?," +
                        "startTime =?," +
                        "endTime=?,"+
                        "fileName=?," +// 11
                        "searchScore=?," +// 12
                        "searchScoreType=?," +// 13
                        "deltaCN=?," +
                        "scan = ? " +
                        "WHERE sequenceCS = ?"
        );
        getScoreStatement = con.prepareStatement(
                "SELECT searchScore, peptideID FROM "+ PEPTIDE_TABLE_NAME +
                        " WHERE sequenceCS = ?"
        );

        getProteinIDFromLocusStatement = con.prepareStatement(
                "SELECT proteinID FROM "+PROTEIN_TABLE_NAME
                +" WHERE Accession = ?;"
        );
        getProteinCountStatement = con.prepareStatement(
                "SELECT count(*) FROM "+PROTEIN_TABLE_NAME+";"
        );

        getSequenceCSIDScoreStatement = con.prepareStatement(
                "SELECT peptideID, sequenceCS, searchScore FROM "+ PEPTIDE_TABLE_NAME+";"
        );
        updateSpectraStatement = con.prepareStatement(
                "UPDATE "+ SPECTRA_TABLE_NAME +
                        " SET peakMZ=?," +
                        "peakIntensity =?," +
                        "massKey = ?" +
                        "WHERE peptideID = ?;"
        );
        addProteinEntryStatement =  con.prepareStatement(
                "INSERT INTO "+PROTEIN_TABLE_NAME+
                        " VALUES(?,?,?); "
        );

        addProteinToPeptideStatement  = con.prepareStatement(
                "INSERT INTO "+PEPTIDE_PROTEIN_INDEX_TABLE +
                        " VALUES(?,?);"
        );


        incrementCopies = con.prepareStatement(
                "UPDATE " + PEPTIDE_TABLE_NAME +
                        " SET copies=copies+1" +
                        " WHERE peptideID = ?;"
        );
//        updateEntry.executeUpdate();


        getSpectraStatement = con.prepareStatement(
                "SELECT * "
                        + "FROM " + SPECTRA_TABLE_NAME + " "
                        + "WHERE peptideID BETWEEN ? AND ?;"
        );
        getSpectraMassStatement = con.prepareStatement(
                "SELECT * "
                        + "FROM " + SPECTRA_TABLE_NAME + " "
                        + "WHERE massKey BETWEEN ? AND ?" +
                        " ORDER BY peptideID;")
        ;
        getProteinIDStatement = con.prepareStatement(
                "SELECT proteinID FROM "+PEPTIDE_PROTEIN_INDEX_TABLE+
                        " WHERE peptideID = ?;"
        );
        getPeptideCountStatement = con.prepareStatement("SELECT count(*) FROM "+
                PEPTIDE_TABLE_NAME +";"
        );
        getPeptideProteinTable = con.prepareStatement("SELECT * FROM "+PEPTIDE_PROTEIN_INDEX_TABLE);

        retreiveSequenceKeys = con.prepareStatement("SELECT sequenceCS FROM "+PEPTIDE_TABLE_NAME);

        retrievePeptideIDs = con.prepareStatement("SELECT peptideID FROM "+PEPTIDE_TABLE_NAME);

        System.out.println(">>>> finished initializing sqlite ");


    }

    private long queryProteinCount() throws SQLException {
        ResultSet rs = getProteinCountStatement.executeQuery();
        if(rs.next())
        {
            return rs.getLong(1);
        }
        return 0;
    }

    private void clearLocalMemory()
    {
        seqKeyIDMap = new HashMap<>();
        seqKeyScoreMap = new HashMap<>();
        seqKeyProteinSetMap = new HashMap<>();
        locusProteinIDMap = new HashMap<>();
    }

    private float queryScore(String key, long id []) throws SQLException {
        getScoreStatement.setString(1,key);
        ResultSet rs = getScoreStatement.executeQuery();
        if(rs.next())
        {
            id[0] = rs.getLong(2);
            return rs.getFloat(1);
        }
        id[0] = -1;
        return Float.MIN_VALUE;
    }


    private int queryPeptideCount() throws SQLException {
        ResultSet rs = getPeptideCountStatement.executeQuery();
        if(rs.next()) return rs.getInt(1);
        return 0;
    }

    private long queryProteinID(String locus) throws SQLException {
        getProteinIDFromLocusStatement.setString(1,locus);
        ResultSet rs = getProteinIDFromLocusStatement.executeQuery();
        if(rs.next())
        {
            return rs.getInt(1);
        }
        return -1L;
    }

    private long uploadProtein(String locus, String description) throws SQLException {
        addProteinEntryStatement.setLong(1,proteinID);
        addProteinEntryStatement.setString(2,locus);
        addProteinEntryStatement.setString(3,description);
        addProteinEntryStatement.executeUpdate();
        return proteinID++;
    }

    private void insertProteinToPeptideIndex(long peptideID, long proteinID) throws SQLException {
        addProteinToPeptideStatement.setLong(1,peptideID);
        addProteinToPeptideStatement.setLong(2,proteinID);
        addProteinToPeptideStatement.executeUpdate();
    }




    public static String getPeptideTableName() {
        return PEPTIDE_TABLE_NAME;
    }

    public long insertEntry(String seq, float precursorMz, int chargeState, int numPeaks, float retTime, float startTime,
                           float endTime, String fileID, float searchScore, int scoreType, float deltaCn,
                           int scan, String sequenceCS)
            throws SQLException
    {
        //String sequenceCS = seq+chargeState;
        int mzKey = (int)(precursorMz*1000);
        addEntryStatement.setLong(1, peptideID); //peptideID 1
        addEntryStatement.setString(2,seq); //peptideSeq 2
        addEntryStatement.setInt(3,mzKey); //mzkey 3
        addEntryStatement.setFloat(4,precursorMz); //mz 4
        addEntryStatement.setInt(5,chargeState); //cs 5
        addEntryStatement.setInt(6,0); //copies 6
        addEntryStatement.setInt(7,numPeaks); //numPeaks 6
        addEntryStatement.setFloat(8,retTime);  //ret time 7
        addEntryStatement.setFloat(9,startTime); //start time 8
        addEntryStatement.setFloat(10,endTime); //end time 9
        addEntryStatement.setString(11,fileID);//filename 10
        addEntryStatement.setFloat(12,searchScore);//search scorea 11
        addEntryStatement.setInt(13,scoreType);// search score type 12
        addEntryStatement.setFloat(14,deltaCn);// delta cn 13
        addEntryStatement.setString(15,sequenceCS);// sequence cs 14
        addEntryStatement.setInt(16,scan);
        addEntryStatement.executeUpdate();
        return peptideID++;
    }


    public static byte[] convertToByteArr(List<Float> floatList) throws IOException {
        ByteArrayOutputStream bous = new ByteArrayOutputStream();
        DataOutputStream dout = new DataOutputStream(bous);
        for(float d: floatList)
        {
            dout.writeFloat(d);
            //bous.write(Float.floatToIntBits(d));
        }
        dout.close();
        return bous.toByteArray();
    }



    public void insertSpectra(long id, float mass, List<Float> mzList, List<Float> intensityList) throws IOException, SQLException {
        List<Long> mzBinList = new ArrayList<>();
        int massKey = (int)(mass*MZ_KEY_SHIFT);
        byte[] mzByteArr = convertToByteArr(mzList);
        byte[] intByteArr = convertToByteArr(intensityList);

        addSpectraStatement.setLong(1,id);
        addSpectraStatement.setBytes(2,mzByteArr);
        addSpectraStatement.setBytes(3,intByteArr);
        addSpectraStatement.setInt(4,massKey);

        addSpectraStatement.executeUpdate();
    }
    private void incrementCopies(int id) throws SQLException {
        incrementCopies.setInt(1,id);
        incrementCopies.executeUpdate();
    }

    public void querySpectra(int id) throws SQLException, IOException {
        getSpectraStatement.setInt(1,id);
        getSpectraStatement.setInt(2,id);
        final ResultSet rs = getSpectraStatement.executeQuery();
        if(rs.next())
        {
            InputStream is = rs.getBinaryStream(2);
            byte [] buffer = new byte[Float.BYTES];
            ByteArrayInputStream bins ;

            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                System.out.print(f+" ");
            }
        }
    }

    public static String removePrePostFix(String sequence)
    {
        String preResult = sequence.substring(0,sequence.lastIndexOf("."));
        return preResult.substring(preResult.indexOf(".")+1,preResult.length());
    }



    public TLongHashSet queryProteinSet(long peptideID) throws SQLException {
        TLongHashSet result = new TLongHashSet();
        getProteinIDStatement.setLong(1,peptideID);
        ResultSet rs = getProteinIDStatement.executeQuery();
        while(rs.next())
        {
            result.add(rs.getLong(1));
        }
        return result;
    }

    public void querySpectra(int mass1, int mass2,   List<LibrarySpectra> spectraList) throws SQLException, IOException {
        int massKey1 = mass1;
        int massKey2 = mass2;
        getSpectraMassStatement.setInt(1,massKey1);
        getSpectraMassStatement.setInt(2,massKey2);
        final ResultSet rs = getSpectraMassStatement.executeQuery();
        int i=0;
        while(rs.next())
        {
            int id = rs.getInt(1);

            InputStream is = rs.getBinaryStream(2);
            byte [] buffer = new byte[Float.BYTES];
            ByteArrayInputStream bins ;
            TFloatArrayList mzList = new TFloatArrayList();
            TFloatArrayList intList = new TFloatArrayList();

            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                mzList.add(f);
            }
            is = rs.getBinaryStream(3);
            buffer = new byte[Float.BYTES];
            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                intList.add(f);
            }
            int mass = rs.getInt(4);
            if(id!=spectraList.get(i).id)
            {
                System.out.println(id+ " "+spectraList.get(i).id+" "+mass);
                System.out.println("Mismatch !");
                System.exit(-1);
            }

            spectraList.get(i).setMzList(mzList);
            spectraList.get(i).setIntensityList(intList);
            i++;
        }
    }

    public void querySpectra(float mass1, float mass2,   List<LibrarySpectra> spectraList) throws SQLException, IOException {
        int massKey1 = (int)(mass1*MZ_KEY_SHIFT);
        int massKey2 = (int)(mass2*MZ_KEY_SHIFT);
        getSpectraMassStatement.setInt(1,massKey1);
        getSpectraMassStatement.setInt(2,massKey2);
        final ResultSet rs = getSpectraMassStatement.executeQuery();
        int i=0;
        while(rs.next())
        {
            int id = rs.getInt(1);

            InputStream is = rs.getBinaryStream(2);
            byte [] buffer = new byte[Float.BYTES];
            ByteArrayInputStream bins ;
            TFloatArrayList mzList = new TFloatArrayList();
            TFloatArrayList intList = new TFloatArrayList();

            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                mzList.add(f);
            }
            is = rs.getBinaryStream(3);
            buffer = new byte[Float.BYTES];
            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                intList.add(f);
            }
            if(id!=spectraList.get(i).id)
            {
                System.out.println(id+ " "+spectraList.get(i).id);
                System.out.println("Mismatch !");
                System.exit(-1);
            }

            spectraList.get(i).setMzList(mzList);
            spectraList.get(i).setIntensityList(intList);
            i++;
        }
    }




    public void querySpectra(int id, List<Float> mzList, List<Float>  intList) throws SQLException, IOException {
        getSpectraStatement.setInt(1,id);
        getSpectraStatement.setInt(2,id);
        final ResultSet rs = getSpectraStatement.executeQuery();
        while(rs.next())
        {
            InputStream is = rs.getBinaryStream(2);
            byte [] buffer = new byte[Float.BYTES];
            ByteArrayInputStream bins ;

            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                mzList.add(f);
            }
            is = rs.getBinaryStream(3);
            buffer = new byte[Float.BYTES];
            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                intList.add(f);
            }
        }
    }

    public float queryEntryGetScore(String sequenceCS) throws SQLException {
        getEntrySequenceCSStatement.setString(1,sequenceCS);
        final ResultSet rs = getEntrySequenceCSStatement.executeQuery();
        if(rs.next())
        {
            return rs.getFloat(12);
        }
        return -1;
    }

    private void updateEntry(String sequenceCS, float precursorMZ, int numPeaks,int copies, float retTime,
                             float startTime, float endTime, String fileName, float searchScore,
                             int scoreType, float deltaCN, int scan) throws SQLException {
        int mzKey = (int)(precursorMZ*MZ_KEY_SHIFT);
        updateEntry.setInt(1,mzKey);
        updateEntry.setFloat(2,precursorMZ);
        updateEntry.setInt(3,numPeaks);
      //  updateEntry.setInt(4,copies);

        updateEntry.setFloat(4,retTime);
        updateEntry.setFloat(5,startTime);
        updateEntry.setFloat(6,endTime);
        updateEntry.setString(7,fileName);
        updateEntry.setFloat(8,searchScore);
        updateEntry.setInt(9,scoreType);
        updateEntry.setFloat(10,deltaCN);
        updateEntry.setString(12,sequenceCS);
        updateEntry.setInt(11,scan);
        updateEntry.executeUpdate();
    }
    private void updateSpectra(long id, List<Float> mzList, List<Float>  intList,float mass) throws IOException, SQLException {
        List<Long> mzBinList = new ArrayList<>();

        byte[] mzByteArr = convertToByteArr(mzList);
        byte[] intByteArr = convertToByteArr(intList);
        int massKey = (int)(mass*MZ_KEY_SHIFT);
        updateSpectraStatement.setLong(4,id);
        updateSpectraStatement.setInt(3,massKey);
        updateSpectraStatement.setBytes(1,mzByteArr);
        updateSpectraStatement.setBytes(2,intByteArr);
        updateSpectraStatement.executeUpdate();
    }

    public List<LibrarySpectra> queryEntry(float start, float end) throws SQLException, IOException {
        List<LibrarySpectra> libList = new ArrayList<>();
        queryEntry(start,end,libList);
        querySpectra(start,end);
        return libList;
    }


    public int[] queryEntry(float start, float end, List<LibrarySpectra> spectraList) throws SQLException {
        int key1 = (int)(start*MZ_KEY_SHIFT);
        int key2 = (int)(end*MZ_KEY_SHIFT);
        int idarr [] = new int[2];
        return queryEntry(key1,key2,spectraList);
    }


    public int[] queryEntry(int start, int end, List<LibrarySpectra> spectraList) throws SQLException {
        int key1 = start;
        int key2 = end;
        int idarr [] = new int[2];
        getEntryByMassRangeStatement.setInt(1,key1);
        getEntryByMassRangeStatement.setInt(2,key2);
        final ResultSet rs = getEntryByMassRangeStatement.executeQuery();
        int i=0;
        long prevId =-1;
        LibrarySpectra prevsSpectra = null;
        while (rs.next()) {
            if(i==0)idarr[0] =rs.getInt(1);
            else idarr[1] = rs.getInt(1);
            long id= rs.getLong(1);
            String accession = rs.getString(20);
            String proteinDescription = rs.getString(21);
           /* if(id==6732)
            {
                System.out.println(">>>");
            }*/

            if(id!=prevId)
            {
                String seq = rs.getString(2);

                float xcorr = rs.getFloat(12);
                float mass = rs.getFloat(4);
                int cs = rs.getInt(5);
                float retTime = rs.getFloat(8);
                float deltaCN = rs.getFloat(14);
                String key = rs.getString(15);
                String fileName = rs.getString(11);
                int scan = rs.getInt(16);

                LibrarySpectra librarySpectra = new LibrarySpectra(seq,cs,mass,retTime, xcorr, deltaCN,key,fileName,scan,id, accession, proteinDescription);
                spectraList.add(librarySpectra);
                prevsSpectra =librarySpectra;
            }
            else prevsSpectra.addProtein(accession,proteinDescription);
            prevId=id;


            i++;
        }
        rs.close();
        return idarr;
    }
    public void readDTASelect(String dtaSelectFile) throws JDOMException, SQLException, IOException {
        String path = dtaSelectFile.substring(0,dtaSelectFile.lastIndexOf(File.separatorChar));
        readDTASelect(dtaSelectFile,path);
    }

    public void readDTASelectCreateCSV(String dtaSelectFile, String workingDirectory,
                                       String outPeptideCSV, String outProteinCSV, Map<String,Long> seqKeyMap)
            throws IOException {
        boolean start = false;
        indexedMap = IndexUtil.createIndexedHt(workingDirectory,"ms2");

        // initializeLocalMaps();
        int peptideRowLength=0;
        int proteinRowLength =0;
        BufferedReader br = new BufferedReader(new FileReader(dtaSelectFile));
        BufferedWriter peptideBW = new BufferedWriter(new FileWriter(outPeptideCSV,true));
        String line;
        String locus = "";
        String description = "";
        String accession = " ";


        //peptideBW.append(CSV_HEADER);
        peptideBW.newLine();
        while((line=br.readLine())!=null)
        {
            if(!start && line.startsWith("Locus"))
            {
                String arr [] = line.split("\t");
                proteinRowLength = arr.length;
                processProteinHeader(arr);
            }
            else if(!start && line.startsWith("Unique") )
            {
                String arr [] = line.split("\t");
                peptideRowLength =arr.length;
                start = true;
                processPeptideHeader(arr);
                continue;
            }
            if(start)
            {
                String arr [] = line.split("\t");
                if(arr.length==1)
                {
                    start = false;
                    continue;
                }
                if(arr[1].startsWith("Proteins"))
                {
                    start = false;
                    continue;
                }
                if(arr.length ==peptideRowLength)
                {
                    processPeptideLineCSV(arr,peptideBW,accession,description,seqKeyMap);
                    peptideBW.newLine();
                    peptideBW.flush();
                }//
                else if(arr.length<=proteinRowLength)
                {
                    locus = arr[locusIndex];
                    description = arr[descriptiveIndex];
                    accession  = getSequestLikeAccession(locus);
                    // System.out.println(getSequestLikeAccession(locus));
                }

            }
        }
        peptideBW.close();
        br.close();
    }



    public void readDTASelect(String dtaSelectFile,String path) throws IOException, JDOMException, SQLException {
        boolean start = false;
         indexedMap = IndexUtil.createIndexedHt(path,"ms2");
         seqKeyScoreMap = new HashMap<>();
         seqKeyIDMap = new HashMap<>();
         peptideID = queryPeptideCount();
        proteinID = queryProteinCount();
       // initializeLocalMaps();
         int peptideRowLength=0;
        int proteinRowLength =0;
        BufferedReader br = new BufferedReader(new FileReader(dtaSelectFile));
        String line;
        String locus = "";
        String description = "";
        String accession = " ";

        while((line=br.readLine())!=null)
        {
            if(!start && line.startsWith("Locus"))
            {
                String arr [] = line.split("\t");
                proteinRowLength = arr.length;
                processProteinHeader(arr);
            }
            else if(!start && line.startsWith("Unique") )
            {
                String arr [] = line.split("\t");
                peptideRowLength =arr.length;
                start = true;
                processPeptideHeader(arr);
                continue;
            }
            if(start)
            {
                String arr [] = line.split("\t");
                if(arr.length==1)
                {
                    start = false;
                    continue;
                }
                if(arr[1].startsWith("Proteins"))
                {
                    start = false;
                    continue;
                }
                if(arr.length ==peptideRowLength)
                {
                       processPeptideLine(arr,path,accession,description);
                }
                else if(arr.length<=proteinRowLength)
                {
                   locus = arr[locusIndex];
                   description = arr[descriptiveIndex];
                    accession  = getSequestLikeAccession(locus);
                   // System.out.println(getSequestLikeAccession(locus));
                }

            }
        }
        br.close();
    }
    private long readAndUpdateSpectra(String seqKey, String filename, String path, float score,
                                      float mass, float deltacn, int scan, long id) throws SQLException, IOException {
        //seqKeyScoreMap.put(seqKey,score);
        //String fpath = path+filename;
        IndexedFile ifile = indexedMap.get(filename);
        List<Float> mzList = new ArrayList<>();
        List<Float> intList = new ArrayList<>();
        //long id = seqKeyIDMap.get(seqKey);
        float retTime = readSpectraFromMS2(ifile,scan,mzList,intList);
        updateEntry(seqKey,mass,mzList.size(),0,retTime,0,0,filename,
                score,0,deltacn,scan);
        updateSpectra(id,mzList,intList,mass);
        return id;
    }

    private void processPeptideLineCSV(String[] peptideArr,Appendable ap, String accession, String description, Map<String, Long> seqKeyMap) throws IOException {
        float xcorr = Float.parseFloat(peptideArr[xcorrIndex]);
        String sequence = peptideArr[seqIndex];
        float mass = Float.parseFloat(peptideArr[theoreticalMassIndex-1]);
        String fileString = peptideArr[fileNameIndex];
        String regex = "\\.";
        String [] fileArr = fileString.split(regex);
        String fileName=  fileArr[0]+".ms2";
        int cs = Integer.parseInt(fileArr[3]);
        int scan = Integer.parseInt(fileArr[2]);
        long peptideID=-1;
        String seqKey = removePrePostFix(sequence)+cs;

       /* if(seqKeyIDMap.containsKey(seqKey))
        {
            peptideID = seqKeyIDMap.get(seqKey);
        }
        else
        {
            peptideID = this.peptideID++;
            seqKeyIDMap.put(seqKey,peptideID);
        }*/
        if(seqKeyMap.containsKey(seqKey))
        {
            peptideID = seqKeyMap.get(seqKey);
        }
        else
        {
            peptideID = LibraryCSV.peptideID++;
            seqKeyMap.put(seqKey, peptideID);
        }


        String fpath = fileName;
        IndexedFile ifile = indexedMap.get(fpath);
        List<Float> mzList = new ArrayList<>();
        List<Float> intList = new ArrayList<>();
        float retTime = readSpectraFromMS2(ifile,scan,mzList,intList);


        StringBuilder sb = new StringBuilder();
        sb.append(peptideID).append("\t");
        sb.append(fpath).append("\t");
        sb.append(sequence).append("\t");
        sb.append(cs).append("\t");
        sb.append(retTime).append("\t");
        sb.append(mass).append("\t");
        sb.append(mzList.size()).append("\t");
        for(float f: mzList)
        {
            sb.append(f).append(",");
        }
        sb.append("\t");
        for(float f: intList)
        {
            sb.append(f).append(",");
        }
        sb.append("\t");
        ap.append(sb);



    }



    private void processPeptideLine(String[] peptideArr,String path, String accession, String description) throws IOException, SQLException {
        float xcorr = Float.parseFloat(peptideArr[xcorrIndex]);
        float deltaCN = Float.parseFloat(peptideArr[deltaCNIndex]);
        String sequence = peptideArr[seqIndex];
        float theoreticalMass = Float.parseFloat(peptideArr[theoreticalMassIndex]);
        float mass = Float.parseFloat(peptideArr[theoreticalMassIndex-1]);
        String fileString = peptideArr[fileNameIndex];
        String regex = "\\.";
        String [] fileArr = fileString.split(regex);
        String fileName=  fileArr[0]+".ms2";
        int cs = Integer.parseInt(fileArr[3]);
        int scan = Integer.parseInt(fileArr[2]);
        long peptideID=-1;
        String seqKey = removePrePostFix(sequence)+cs;


        float libraryScore = Float.MIN_VALUE;
        TLongHashSet proteinSet = null;
        boolean peptideExistsInLibrary = true;

        /*
               Check if entry exists in local memory;
               if not check it it exists in library;
                    if in library load into memory
                    else  set flag for upload
         */

        if(!seqKeyScoreMap.containsKey(seqKey))
        {
            long [] idarr = new  long[1];
             libraryScore = queryScore(seqKey,idarr);

             if((idarr[0]!=-1))
             {
                 peptideID = idarr[0];
                seqKeyScoreMap.put(seqKey,libraryScore);
                seqKeyIDMap.put(seqKey,peptideID);
                proteinSet = queryProteinSet(peptideID);
                seqKeyProteinSetMap.put(seqKey,proteinSet);
             }
             else peptideExistsInLibrary = false;
        }
        else
        {
            peptideID = seqKeyIDMap.get(seqKey);
            libraryScore = seqKeyScoreMap.get(seqKey);
            proteinSet = seqKeyProteinSetMap.get(seqKey);
        }

        if(!peptideExistsInLibrary)
        {
            String fpath = fileName;
            IndexedFile ifile = indexedMap.get(fpath);
            List<Float> mzList = new ArrayList<>();
            List<Float> intList = new ArrayList<>();
            float retTime = readSpectraFromMS2(ifile,scan,mzList,intList);
            peptideID = insertEntry(sequence,mass,cs,mzList.size(),retTime,0,0,fileName,
                    xcorr,0,deltaCN,scan,seqKey);
            seqKeyIDMap.put(seqKey,peptideID);
            proteinSet = new TLongHashSet();
            seqKeyProteinSetMap.put(seqKey, proteinSet);
            seqKeyScoreMap.put(seqKey,xcorr);

            insertSpectra(peptideID,mass,mzList,intList);
        }
        /*
          check if new score is higher than library score
                upload if true
        */
        else if(xcorr>libraryScore) {
            seqKeyScoreMap.put(seqKey,xcorr);
            readAndUpdateSpectra(seqKey, fileName, path, xcorr, mass, deltaCN, scan, peptideID);
        }
        /*
                Check if accession exists in local memory
                    if not check if exists in library
                        if not add protein to protein table
                            flag to add to peptide
                    else
                        add to local memory
         */
        long proteinID = -1;
        boolean addProteinToPeptide = false;
        if(!locusProteinIDMap.containsKey(seqKey))
        {
            proteinID = queryProteinID(accession);
            if(proteinID!=-1)
            {
                locusProteinIDMap.put(accession,proteinID);
            }
            else
            {
                addProteinToPeptide = true;
                proteinID = uploadProtein(accession,description);
            }
        }
        else proteinID = locusProteinIDMap.get(seqKey);
        if(addProteinToPeptide || !proteinSet.contains(proteinID))
        {
            proteinSet.add(proteinID);
            //System.out.println(peptideID+ " " +proteinID);
            insertProteinToPeptideIndex(peptideID,proteinID);
        }



        /*
        //check if entry exists in local memory/
        if(seqKeyScoreMap.containsKey(seqKey))
        {
             libraryScore = seqKeyScoreMap.get(seqKey);
            if(xcorr>libraryScore)
            {
                peptideID = readAndUpdateSpectra(seqKey,fileName,path,xcorr,mass,deltaCN,scan);
            }
        }
        else
        {
            long [] idarr = new  long[1];
             libraryScore = queryScore(seqKey,idarr);

            if(xcorr>libraryScore)
            {
                readAndUpdateSpectra(seqKey,fileName,path,xcorr,mass,deltaCN,scan);
            }
            if(idarr[0]!=-1)
            {
                seqKeyScoreMap.put(seqKey,libraryScore);
                seqKeyIDMap.put(seqKey,idarr[0]);
                //seqKeyProteinSetMap.put(seqKey,)
            }
            else
            {

                seqKeyScoreMap.put(seqKey,xcorr);
                String fpath = fileName;
                IndexedFile ifile = indexedMap.get(fpath);
                List<Float> mzList = new ArrayList<>();
                List<Float> intList = new ArrayList<>();
                float retTime = readSpectraFromMS2(ifile,scan,mzList,intList);
                 peptideID = insertEntry(sequence,mass,cs,mzList.size(),retTime,0,0,fileName,
                        xcorr,0,deltaCN,scan,seqKey);
                seqKeyIDMap.put(seqKey,peptideID);
                insertSpectra(peptideID,mass,mzList,intList);

            }
        }*/

    }

    private float readSpectraFromMS2(IndexedFile ifile, int scan,List<Float> mzList, List<Float> intList)
            throws IOException {
        long pos = ifile.getPositionByScan(scan);
        RandomAccessFile rfile = new RandomAccessFile(ifile.getFileName(),"r");
        rfile.seek(pos);
        float retTime = -1;
        String line=rfile.readLine();

        while((line=rfile.readLine())!=null && !line.startsWith("S")){
            if(line.startsWith("I\tRetTime\t"))
            {
                retTime = Float.parseFloat( line.split("\t")[2]);

            }
            else if(Character.isDigit(line.charAt(0)))
            {
                String[] arr = line.split(" ");
                mzList.add(Float.parseFloat(arr[0]));
                intList.add(Float.parseFloat(arr[1]));
            }
        }
        rfile.close();
        return retTime;
    }

    private void processProteinHeader(String [] proteinHeader)
    {
        int i=0;
        for(String str:proteinHeader)
        {
            if(str.contentEquals( "Descriptive Name"))
            {
                descriptiveIndex =i;
            }
            else if(str.contentEquals( "Locus"))
            {
                locusIndex = i;
            }
            i++;
        }
    }
    private void processPeptideHeader(String [] peptideHeader )
    {
        int i=0;
        for(String str:peptideHeader)
        {
            if(str.contentEquals( "Sequence"))
            {
                seqIndex =i;
            }
            else if(str.contentEquals( "FileName"))
            {
                fileNameIndex = i;
            }
            else if(str.contentEquals("CalcM+H+"))
            {
                theoreticalMassIndex = i;
            }
            else if(str.contentEquals("XCorr"))
            {
                xcorrIndex=i;
            }
            else if(str.contentEquals("DeltCN"))
            {
                deltaCNIndex =i;
            }
            i++;
        }
    }

    public List<LibrarySpectra> querySpectra(float mzStart, float mzEnd) throws SQLException, IOException {
        List<LibrarySpectra> result = new ArrayList<>();
        int idArr[] =queryEntry(mzStart,mzEnd,result);
        querySpectra(mzStart,mzEnd,result);
        return result;
    }

    public List<LibrarySpectra> querySpectra(int mzStart, int mzEnd) throws SQLException, IOException {
        List<LibrarySpectra> result = new ArrayList<>();
        int idArr[] =queryEntry(mzStart,mzEnd,result);
        querySpectra(mzStart,mzEnd,result);
        return result;
    }





    public static void main(String [] args) throws SQLException, IOException, JDOMException {
        //testLibCreation();
        //testReadingDTASelect();
        //testReadLibrary();
//        testReadLib2();
        String libAddress = args[0];
        File f = new File(libAddress);
        if(f.exists()) f.delete();
        LibraryIndexer libraryIndexer = new LibraryIndexer(libAddress);
        for(int i=1; i<args.length; i++)
        {
            System.out.println("reading "+args[i]);
            libraryIndexer.readDTASelect(args[i]);
        }

      /*  System.out.println("Working Directory = " +
                System.getProperty("user.dir"));*/
        //testReadingDTASelect2();

    }

    public static void testReadLib2() throws SQLException, IOException {
        String address="/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/testLibDTASelect.lib";
        LibraryIndexer libraryIndexer = new LibraryIndexer(address);
        List<LibrarySpectra> spectraList = libraryIndexer.querySpectra(1813.0f,1814.f);
        for(LibrarySpectra spectra: spectraList)
        {
            System.out.println(spectra);
            System.out.println(spectra.printSpectra());
        }
    }


    public static void testReadingDTASelect2() throws SQLException, IOException, JDOMException {
        String address="/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/testLibDTASelect.lib";



        LibraryIndexer libraryIndexer = new LibraryIndexer(address);
        /*libraryIndexer.readDTASelect("/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/DTASelect-filter.txt",
                "/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/","");*/
        libraryIndexer.readDTASelect("/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/DTASelect-filter3.txt",
                "/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/");
    }



    public static void testReadingDTASelect() throws SQLException, IOException, JDOMException {
        String address="/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/testLibDTASelect.lib";
        File f = new File(address);
        if(f.exists()) f.delete();


        LibraryIndexer libraryIndexer = new LibraryIndexer(address);
        libraryIndexer.readDTASelect("/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/DTASelect-filter.txt",
                "/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/");
        libraryIndexer.readDTASelect("/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/DTASelect-filter2.txt",
                "/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/");
        libraryIndexer.readDTASelect("/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/DTASelect-filter3.txt",
                "/home/yateslab/project_data/LibraryBuilder/2803sqliteBuild/");
    }

    public void queryPeptideProteinTable(List<Integer> peptideIDList, List<Integer> proteinIDList) throws SQLException {
        ResultSet rs = getPeptideProteinTable.executeQuery();
        while(rs.next())
        {
            peptideIDList.add(rs.getInt(1));
            proteinIDList.add(rs.getInt(2));
        }
    }

    public List<Integer> queryPeptideIDS() throws SQLException {
        List<Integer> peptideIDList = new ArrayList<>();
        ResultSet rs = retrievePeptideIDs.executeQuery();
        while(rs.next())
        {
            peptideIDList.add(rs.getInt(1));
        }
        return peptideIDList;
    }

    public List<String> queryPeptideKeys() throws SQLException {
        List<String> peptideKeys = new ArrayList<>();
        ResultSet rs = retreiveSequenceKeys.executeQuery();
        while(rs.next())
        {
            peptideKeys.add(rs.getString(1));
        }
        return peptideKeys;
    }


}
