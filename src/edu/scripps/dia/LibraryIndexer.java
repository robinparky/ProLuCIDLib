package edu.scripps.dia;

import com.google.common.collect.TreeMultimap;
import edu.scripps.dia.util.IndexUtil;
import edu.scripps.dia.util.IndexedFile;
import edu.scripps.dia.util.SpectrumReader;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;
import gnu.trove.TFloatArrayList;
import gnu.trove.TLongHashSet;
import org.jdom.JDOMException;
import org.sqlite.SQLiteConfig;

import java.io.*;
import java.nio.ByteBuffer;
import java.sql.*;
import java.util.*;

import static edu.scripps.pms.util.seq.Fasta.getSequestLikeAccession;
import static rpark.statistics.Outlier.removeOutlier;

/**
 * Created by Titus Jung titusj@scripps.edu on 3/26/18.
 */
public class LibraryIndexer {

    private PreparedStatement addEntryStatement;
    private PreparedStatement addSpectraStatement;
    private PreparedStatement getEntryByMassRangeStatement;
    private PreparedStatement getEntryByMassRangeGTStatement;
    private PreparedStatement getEntrySequenceCSStatement;
    private PreparedStatement updateEntry;
    private PreparedStatement getSpectraStatement;
    private PreparedStatement incrementCopies;
    private PreparedStatement getSequenceCSIDScoreStatement;
    private PreparedStatement getScoreStatement;
    private PreparedStatement getNextHighestMassIndex;
    private PreparedStatement getEntryByMassIndex;
    protected PreparedStatement updateSpectraStatement;
    protected Connection con;
    private long peptideID =0;
    private long proteinID = 0;
    private long spectraID = 0;
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
    private PreparedStatement getSpectraCountStatement;
    private PreparedStatement insertSpectraMetaTable;
    private PreparedStatement createSpectraMetaTable;
    private PreparedStatement getNPeptidesStatement;
    private PreparedStatement addIsDecoyColumnToSpectraMetaTable;
    private PreparedStatement addHasDecoyColumnToSpectraMetaTable;
    private PreparedStatement fillMetaTable;
    private PreparedStatement getNMetaForwardSpectraStatement;
    private PreparedStatement getNMetaAllSpectraStatement;
    private PreparedStatement getRetentionTimesForPeptide =null;
    private PreparedStatement setRetTimeRangeForPeptideID =null;
    private PreparedStatement setRetentionTimeList = null;
    private PreparedStatement getRetentionTimeList = null;
    private PreparedStatement setRetTimeRangeForSeqKey = null;


    private boolean uploadBestScoringXcorr = true;
    private boolean uploadAllSpectra = false;
    private boolean appendAllSpectra = false;


    private Set<String> fileNameScanNumSet = new HashSet<>();
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
    private Map<String, List<Double>>   seqKeyRetTimeMap = new HashMap<>();
    private Set<String> fileStringSet = new HashSet<>();
    public static final float NO_FLOAT_VALUE = -2000;
    public static final String CSV_HEADER =
            "H\tPeptideID\tFilePath\tSequence\tChargeState\tRetentionTime\tMass\tNumberOfPeaks\tMassList\tIntensityList";

    public final static int MZ_KEY_SHIFT = 1000;
    public static final float DECOY_MASS_SHIFT_FLOAT =8;

    public static final int DECOY_MASS_SHIFT =(int) DECOY_MASS_SHIFT_FLOAT*MZ_KEY_SHIFT;



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
            "sequenceCS," +// 15
            "scan," +// 16
            "isDecoy," +// 17
            "hasDecoy, " +// 18
            "retentionTimeList) ";




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
            "scan INTEGER," +
            "isDecoy INTEGER," +
            "hasDecoy INTEGER," +
            "retentionTimeList TEXT) ";

    public static final String PEPTIDE_TABLE_NAME = "PeptideTable";

    public static final String SPECTRA_TABLE_NAME = "SpectraTable";
    public static final String spectraShema = " (peptideID,peakMZ,peakIntensity,massKey, retTime, fileName, scanNumber)";
    public static final String spectraShemaCreation
            = " (peptideID INTEGER," +
            "peakMZ BLOB," +
            "peakIntensity BLOB," +
            "massKey Integer," +
            "retTime FLOAT," +
            "fileName TEXT," +
            "scanNumber Integer)";

    public static final String SPECTRA_META_TABLE_NAME = "SpectraMetaTable";
    public static final String SPECTRA_META_TABLE_SCHEMA ="(spectraID, massKey)";
    public static final String SPECTRA_META_TABLE_CREATION0 = "(spectraID INTEGER, massKey INTEGER," +
            "chargeState INTEGER, isDecoy INTEGER DEFAULT 0, hasDecoy INTEGER DEFAULT 0, diff FLOAT DEFAULT 0)";



    public static final String SPECTRA_META_TABLE_CREATION =
            "CREATE TABLE "+ SPECTRA_META_TABLE_NAME + " AS SELECT "+PEPTIDE_TABLE_NAME+".ROWID AS spectraID, "+
                    PEPTIDE_TABLE_NAME+".massKey, "+PEPTIDE_TABLE_NAME+".chargeState FROM " +SPECTRA_TABLE_NAME+", "+
                    PEPTIDE_TABLE_NAME + " WHERE " +SPECTRA_TABLE_NAME+".peptideID = "+SPECTRA_TABLE_NAME+".peptideID";

    private PreparedStatement unsertSpectraMetaTableStatement;
    public static final String INSERT_SPECTRA_META_TABLE  = "INSERT INTO "+SPECTRA_META_TABLE_NAME + " VALUES (?,?,?,?,?,?);";

    private PreparedStatement getNumberOfDecoysStatement;
    public static final String COUNT_NUMBER_OF_DECOYS_STATEMENT = "SELECT count(*) FROM "+SPECTRA_META_TABLE_NAME+
            " WHERE isDecoy = 1";
    private PreparedStatement updateHasDecoyStateStatement;
    public static final String UPDATE_HAS_DECOY_STATE_STATEMENT = "UPDATE "+SPECTRA_META_TABLE_NAME + " SET hasDecoy =1"
            + " WHERE massKey>=? AND massKey<=? AND isDecoy = 0 ;";
    private PreparedStatement getSpectraMassWithDecoysStatement;
    public static final String GET_SPECTRA_MASS_WITH_DECOYS_STATEMENT =
            "SELECT  SpectraMetaTable.massKey, SpectraTable.peptideID,SpectraMetaTable.chargeState," +
                    "SpectraMetaTable.isDecoy, SpectraTable.peakMZ, SpectraTable.peakIntensity FROM SpectraMetaTable " +
                    "INNER JOIN SpectraTable ON SpectraTable.rowID = SpectraMetaTable.spectraId " +
                    "WHERE SpectraMetaTable.massKey>? AND SpectraMetaTable.massKey<? ORDER BY SpectraMetaTable.massKey;";
    private PreparedStatement getDecoySpectraStatement;
    public static final String GET_DECOY_SPECTRA_STATEMENT =
            "SELECT  SpectraMetaTable.massKey, SpectraMetaTable.chargeState," +
                    "SpectraTable.peakMZ, SpectraTable.peakIntensity, SpectraMetaTable.diff  FROM SpectraMetaTable " +
                    "INNER JOIN SpectraTable ON SpectraTable.rowID = SpectraMetaTable.spectraId " +
                    "WHERE SpectraMetaTable.massKey>? AND SpectraMetaTable.massKey<?  AND SpectraMetaTable.isDecoy =1 " +
                    "ORDER BY SpectraMetaTable.massKey;";
    private PreparedStatement getDecoySpectraWithPeptideAndProteinStatement = null;
    public static final String GET_DECOY_ENTRY_WITH_PEPTIDE_AND_PROTEIN_STATEMENT =
            "SELECT  SpectraMetaTable.rowID, SpectraMetaTable.massKey, SpectraMetaTable.chargeState," +
                    "peptideSeq, accession, description, PeptideTable.sequenceCS FROM SpectraMetaTable " +
                    "INNER JOIN SpectraTable ON SpectraTable.rowID = SpectraMetaTable.spectraId " +
                    "INNER JOIN PeptideTable ON PeptideTable.peptideID = SpectraTable.peptideId " +
                    "INNER JOIN PeptideProteinIndexTable ON PeptideProteinIndexTable.peptideID = PeptideTable.peptideID " +
                    "INNER JOIN ProteinTable ON ProteinTable.proteinID = PeptideProteinIndexTable.proteinID " +
                    "WHERE SpectraMetaTable.massKey BETWEEN ? and ?   AND SpectraMetaTable.isDecoy =1 " +
                    "ORDER BY SpectraMetaTable.massKey;";
    private PreparedStatement getAllSpectraInformation;
    public static final String GET_ALL_SPECTRA_INFO_STATMENT =
            "SELECT SpectraTable.rowId, SpectraTable.massKey, PeptideTable.chargestate FROM SpectraTable, peptideTable " +
                    "WHERE SpectraTable.peptideId = peptideTable.peptideId and SpectraTable.peptideId >=0 and peptidetable.isdecoy =0;";

    private PreparedStatement clearSpectraMetaTable;
    public static final String CLEAR_SPECTRA_META_TABLE = "DELETE FROM "+SPECTRA_META_TABLE_NAME;

    private Map<Integer, List<Double>> retentionTimeMap = new HashMap<>();

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
        System.out.println(":: db at "+dbAddress);
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
        con.prepareStatement("CREATE TABLE IF NOT EXISTS "+SPECTRA_META_TABLE_NAME
                +" "+SPECTRA_META_TABLE_CREATION0).execute();
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


        getNumberOfDecoysStatement = con.prepareStatement(COUNT_NUMBER_OF_DECOYS_STATEMENT);

        addEntryStatement
                = con.prepareStatement("INSERT INTO " + getPeptideTableName() + peptideSchema
                + " VALUES (?, ?, ?, ?, ?, ?, ?, ?,?, ?,?, ?, ?,?,?,?,?,?,?);") ;

        getEntryByMassRangeStatement = con.prepareStatement(
                "SELECT * "
                        + "FROM " + PEPTIDE_TABLE_NAME + " "
                        + "INNER JOIN "+ PEPTIDE_PROTEIN_INDEX_TABLE
                        + " ON " +PEPTIDE_TABLE_NAME+".peptideID="+PEPTIDE_PROTEIN_INDEX_TABLE+".peptideID "
                        + "INNER JOIN "+PROTEIN_TABLE_NAME
                        +" ON "+PEPTIDE_PROTEIN_INDEX_TABLE+".proteinID ="+PROTEIN_TABLE_NAME+".proteinID "
                        + "WHERE massKey BETWEEN ? AND ? AND isDecoy = 0" +
                        " ORDER BY peptideID;")
        ;

    //   updateHasDecoyStateStatement = con.prepareStatement(UPDATE_HAS_DECOY_STATE_STATEMENT);
        getDecoySpectraStatement = con.prepareStatement(GET_DECOY_SPECTRA_STATEMENT);
        getAllSpectraInformation = con.prepareStatement(GET_ALL_SPECTRA_INFO_STATMENT);
        getEntrySequenceCSStatement= con.prepareStatement(
                "SELECT * "
                        + "FROM " + PEPTIDE_TABLE_NAME + " "
                        + "WHERE sequenceCS = ?;");
        updateEntry = con.prepareStatement(
                "UPDATE "+ PEPTIDE_TABLE_NAME +
                        " SET massKey =?," + //1
                        "precursorMZ =?," +//2
                        "numPeaks =?," + //3
                 //       "copies =?,"+
                        "retentionTime =?," + //4
                        "startTime =?," +  //5vim o
                        "endTime=?,"+ //6
                        "fileName=?," +// 7
                        "searchScore=?," +// 8
                        "searchScoreType=?," +// 9
                        "deltaCN=?," + // 10
                        "scan = ?, " +//11
                        "isDecoy = ?," +//12
                        "hasDecoy = ?"+//13
                        "WHERE sequenceCS = ?"//14
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

        getSpectraCountStatement = con.prepareStatement(
                "SELECT count(*) FROM "+SPECTRA_TABLE_NAME+";"
        );

        insertSpectraMetaTable = con.prepareStatement(INSERT_SPECTRA_META_TABLE);

        getSequenceCSIDScoreStatement = con.prepareStatement(
                "SELECT peptideID, sequenceCS, searchScore FROM "+ PEPTIDE_TABLE_NAME+";"
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
                "SELECT SpectraTable.* "
                        + "FROM " + SPECTRA_TABLE_NAME + " INNER JOIN PeptideTable ON SpectraTable.peptideID = PeptideTable.peptideID  "
                        + "WHERE SpectraTable.massKey BETWEEN ? AND ? AND  isDecoy = 0 " +
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

        getNextHighestMassIndex = con.prepareStatement("SELECT massKey FROM "+PEPTIDE_TABLE_NAME
                +" WHERE massKey>? ORDER BY massKey LIMIT 1");

        getNMetaForwardSpectraStatement = con.prepareStatement("SELECT rowID,spectraID,massKey,chargeState FROM "
                +SPECTRA_META_TABLE_NAME+" WHERE massKey >= ? AND isDecoy=0 AND hasDecoy =0 " +
                "ORDER BY massKey LIMIT ? ");

        getNMetaAllSpectraStatement  = con.prepareStatement("SELECT rowID,spectraID,massKey,chargeState FROM "
                +SPECTRA_META_TABLE_NAME+" WHERE massKey >= ? AND isDecoy=0 " +
                "ORDER BY massKey LIMIT ? ");

        getSpectraMassWithDecoysStatement = con.prepareStatement(GET_SPECTRA_MASS_WITH_DECOYS_STATEMENT);

        //        con.commit();
        clearSpectraMetaTable = con.prepareStatement(CLEAR_SPECTRA_META_TABLE);

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

    private long querySpectraCount() throws SQLException {
        ResultSet rs = getSpectraCountStatement.executeQuery();
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

    private long insertProtein(String locus, String description) throws SQLException {
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
                            int scan, String sequenceCS) throws SQLException {
        return insertEntry(seq,precursorMz,chargeState,numPeaks,retTime,startTime,endTime,fileID,
                searchScore,scoreType,deltaCn,scan,sequenceCS,false);
    }


    public long insertEntry(String seq, float precursorMz, int chargeState, int numPeaks, float retTime, float startTime,
                           float endTime, String fileID, float searchScore, int scoreType, float deltaCn,
                           int scan, String sequenceCS, boolean isDecoy)
            throws SQLException
    {
        //String sequenceCS = seq+chargeState;
        int decoyCode = isDecoy? 1: 0;

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
        addEntryStatement.setInt(17,decoyCode);
        addEntryStatement.setInt(18,0);
        addEntryStatement.setString(19,"");

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

    public void insertSpectra(long id, float mass, List<Float> mzList, List<Float> intensityList, float retTime,
                              String fileName, int scanNumber) throws IOException, SQLException {
        insertSpectra(id,mass,mzList,intensityList,retTime,fileName,scanNumber,false);
    }

    public void insertSpectra(long id, float mass, List<Float> mzList, List<Float> intensityList, float retTime,
                              String fileName, int scanNumber, boolean isDecoy) throws IOException, SQLException {
        List<Long> mzBinList = new ArrayList<>();
        int massKey = (int)(mass*MZ_KEY_SHIFT);
        byte[] mzByteArr = convertToByteArr(mzList);
        byte[] intByteArr = convertToByteArr(intensityList);
        if(addSpectraStatement == null)
        {
            addSpectraStatement = con.prepareStatement("INSERT INTO " + SPECTRA_TABLE_NAME + " "+spectraShema
                    + "VALUES (?, ?, ?,?,?,?,?);") ;
        }

        addSpectraStatement.setLong(1,id);
        addSpectraStatement.setBytes(2,mzByteArr);
        addSpectraStatement.setBytes(3,intByteArr);
        addSpectraStatement.setInt(4,massKey);
        addSpectraStatement.setFloat(5,retTime);
        addSpectraStatement.setString(6,fileName);
        addSpectraStatement.setInt(7,scanNumber);
       // addSpectraStatement.setBoolean(8, isDecoy);
        //addSpectraStatement.setInt(8,1);
        spectraID++;
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

    public List<LibrarySpectra> queryDecoyEntryWithPeptideProtein(int mass1, int mass2)
            throws SQLException, IOException {
        int massKey1 = mass1;
        int massKey2 = mass2;
        if(getDecoySpectraWithPeptideAndProteinStatement == null)
        {
            String temp = "SELECT  SpectraMetaTable.rowID, SpectraMetaTable.massKey, SpectraMetaTable.chargeState," +
                    "peptideSeq, accession, description, PeptideTable.sequenceCS, PeptideTable.retentionTime, " +
                    "PeptideTable.startTime, PeptideTable.endTime FROM SpectraMetaTable " +
                    "INNER JOIN SpectraTable ON SpectraTable.rowID = SpectraMetaTable.spectraId " +
                    "INNER JOIN PeptideTable ON PeptideTable.peptideID = SpectraTable.peptideId " +
                    "INNER JOIN PeptideProteinIndexTable ON PeptideProteinIndexTable.peptideID = PeptideTable.peptideID " +
                    "INNER JOIN ProteinTable ON ProteinTable.proteinID = PeptideProteinIndexTable.proteinID " +
                    "WHERE SpectraMetaTable.massKey BETWEEN ? and ?   AND SpectraMetaTable.isDecoy =1 " +
                    "ORDER BY SpectraMetaTable.massKey;";
            getDecoySpectraWithPeptideAndProteinStatement = con.prepareStatement(temp);
        }
        getDecoySpectraWithPeptideAndProteinStatement.setInt(1,massKey1);
        getDecoySpectraWithPeptideAndProteinStatement.setInt(2,massKey2);
        List<LibrarySpectra> spectraList = new ArrayList<>();
        final ResultSet rs = getDecoySpectraWithPeptideAndProteinStatement.executeQuery();
        int i=0;
        int prevId = -1;
        LibrarySpectra spectra = null;
        while(rs.next())
        {
            int id = rs.getInt(1);
            if(id == prevId)
            {
                String accession = rs.getString(5);
                String description = rs.getString(6);
                spectra.addProtein(accession,description);

            }
            else
            {
                int massId = rs.getInt(2);
                float mass = massId/MZ_KEY_SHIFT;
                int cs = rs.getInt(3);
                String sequence = rs.getString(4);

                String accession = rs.getString(5);
                String description = rs.getString(6);
                String key = rs.getString(7);
                float retTime = rs.getFloat( 8);
                float startTime = rs.getFloat(9);
                float endTime = rs.getFloat(10);
                spectra = new LibrarySpectra(massId,cs,sequence,  true,key,retTime );
                spectra.setEndTime(endTime);
                spectra.setStartTime(startTime);
                spectra.addProtein(accession,description);
                spectraList.add(spectra);
                prevId = id;
            }
            i++;
        }
        return spectraList;
    }

    private PreparedStatement getUnIdentifiedSpectraStatement = null;

    public List<LibrarySpectra> queryUnIdentifiedSpectra(int mass1, int mass2, List<LibrarySpectra> decoySpectraList) throws SQLException, IOException {
        if(getUnIdentifiedSpectraStatement== null)
        {
            getUnIdentifiedSpectraStatement = con.prepareStatement("SELECT PeptideTable.massKey,  " +
                    "chargeState, startTime, endTime," +
                    "PeptideTable.retentionTime, peakMZ,  peakIntensity, sequenceCS," +
                    "PeptideTable.filename, PeptideTable.scan " +
                    " FROM PeptideTable INNER JOIN SpectraTable ON PeptideTable.peptideID = SpectraTable.peptideID " +
                    "WHERE PeptideTable.massKey>? AND PeptideTable.massKey<?  AND PeptideTable.isDecoy = 1;" );
        }
        getUnIdentifiedSpectraStatement.setInt(1, mass1);
        getUnIdentifiedSpectraStatement.setInt(2, mass2);
        ResultSet rs = getUnIdentifiedSpectraStatement.executeQuery();
        List<LibrarySpectra> spectraList = new ArrayList<>();
        int count = 0;
        while(rs.next())
        {
            int massKey = rs.getInt(1);
            int cs  = rs.getInt(2);
            double startTime = rs.getDouble(3);
            double endTime = rs.getDouble(4);
            float retTime = rs.getFloat(5);
            String key = rs.getString(8);

            String fileName = rs.getString(9);
            String scan = rs.getString(10);

            InputStream is = rs.getBinaryStream(6);
            byte [] buffer = new byte[Float.BYTES];
            ByteArrayInputStream bins ;
            TFloatArrayList mzList = new TFloatArrayList();
            TFloatArrayList intList = new TFloatArrayList();

            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                mzList.add(f);
            }
            is = rs.getBinaryStream(7);
            buffer = new byte[Float.BYTES];
            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                intList.add(f);
            }
            count = count >= decoySpectraList.size()-1 ? 0 : count+1;
            String seq  ;
            String description;

            if(decoySpectraList.size()>0)
            {
                LibrarySpectra decoySpectra = decoySpectraList.get(count);
                 seq  = decoySpectra.sequence;
                 description = fileName+scan;
            }
            else
            {
                String tempSeq = fileName.replaceAll("_","").replaceAll("\\.ms2","")
                        .replaceAll("[0-9]","").toUpperCase();
                seq = "N."+ tempSeq + ".F";
                description = fileName;
            }


            String accession = "_from_"+description;

            LibrarySpectra spectra = new LibrarySpectra(massKey,cs,seq,  true,key,retTime );
            spectra.addProtein(accession,description);
            spectra.setIntensityList(intList);
            spectra.setMzList(mzList);
            spectraList.add(spectra);

        }
        return spectraList;
    }

    public List<LibrarySpectra> queryDecoySpectra(int mass1, int mass2) throws SQLException, IOException {
        int massKey1 = mass1;
        int massKey2 = mass2;
        getDecoySpectraStatement.setInt(1,massKey1);
        getDecoySpectraStatement.setInt(2,massKey2);
        List<LibrarySpectra> spectraList = queryDecoyEntryWithPeptideProtein(mass1,mass2);
        final ResultSet rs = getDecoySpectraStatement.executeQuery();
        int i=0;
        while(rs.next())
        {
            int id = rs.getInt(1);
            float mass = id/MZ_KEY_SHIFT;
            int cs = rs.getInt(2);
            float diff = rs.getFloat(5);
            InputStream is = rs.getBinaryStream(3);
            byte [] buffer = new byte[Float.BYTES];
            ByteArrayInputStream bins ;
            TFloatArrayList mzList = new TFloatArrayList();
            TFloatArrayList intList = new TFloatArrayList();

            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                f+=diff;
                mzList.add(f);
            }
            is = rs.getBinaryStream(4);
            buffer = new byte[Float.BYTES];
            while((is.read(buffer))>0)
            {
                float f =ByteBuffer.wrap(buffer).getFloat();
                intList.add(f);
            }

           // LibrarySpectra spectra = new LibrarySpectra(id,cs, true);
            LibrarySpectra spectra = spectraList.get(i);

            if(id!=spectraList.get(i).massKey)
            {
                System.out.println(id+ " "+spectraList.get(i).id+" "+mass);
                System.out.println("Mismatch At Decoy level !");
                System.exit(-1);
            }

            spectra.setMzList(mzList);
            spectra.setIntensityList(intList);
          //  spectraList.add(spectra);
            i++;
        }
        return spectraList;
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
        updateEntry.setInt(13,0);
        updateEntry.setInt(12,0);

        updateEntry.setString(14,sequenceCS);
        updateEntry.setInt(11,scan);
        updateEntry.executeUpdate();
    }

    private void updateSpectra(long id, List<Float> mzList, List<Float>  intList, float mass, float retTime,
                               String filename, int scanNumber) throws IOException, SQLException {
        List<Long> mzBinList = new ArrayList<>();

        byte[] mzByteArr = convertToByteArr( mzList);
        byte[] intByteArr = convertToByteArr(intList);
        int massKey = (int)(mass*MZ_KEY_SHIFT);

        updateSpectraStatement = con.prepareStatement(
                "UPDATE "+ SPECTRA_TABLE_NAME +
                        " SET peakMZ=?," +
                        "peakIntensity =?," +
                        "massKey = ?," +
                        "retTime = ? ," +
                        "fileName = ?," +
                        "scanNumber = ? "+
                        "WHERE peptideID = ?;"
        );

        updateSpectraStatement.setLong(7,id);
        updateSpectraStatement.setFloat(6,scanNumber);
        updateSpectraStatement.setString(5,filename);

        updateSpectraStatement.setFloat(4,retTime);
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
            String accession = rs.getString(23);
            String proteinDescription = rs.getString(24);
           /* if(id==6732)
            {
                System.out.println(">>>");
            }*/

            if(id!=prevId)
            {
                String seq = rs.getString(2);
                int massKey = rs.getInt(3);
                float xcorr = rs.getFloat(12);
                float mass = rs.getFloat(4);
                int cs = rs.getInt(5);
                float retTime = rs.getFloat(8);
                float startTime = rs.getFloat(9);
                float endTime = rs.getFloat(10);
                float deltaCN = rs.getFloat(14);
                String key = rs.getString(15);
                String fileName = rs.getString(11);
                int scan = rs.getInt(16);

                LibrarySpectra librarySpectra = new LibrarySpectra(seq,cs,mass,retTime, xcorr, deltaCN,key,fileName,
                        scan,id, accession, proteinDescription, massKey);
                librarySpectra.setStartTime(startTime);
                librarySpectra.setEndTime(endTime);
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

    public List<Double> getAllRetentionTimesForPeptide(int id) throws SQLException {
        List<Double> retTimeList =new ArrayList<>();
        if(getRetentionTimesForPeptide == null)
        {
            getRetentionTimesForPeptide = con.prepareStatement("SELECT retTime FROM " +SPECTRA_TABLE_NAME+
                    " where peptideID = ?");
        }
        getRetentionTimesForPeptide.setInt(1,id);
        final ResultSet rs = getRetentionTimesForPeptide.executeQuery();
        while(rs.next())
        {
            double retTime = rs.getDouble(1);
            retTimeList.add(retTime);
        }
        return retTimeList;
    }

    private void setRetTimeRangeForPeptide(int id, double startTime, double endTime) throws SQLException {
        if(setRetTimeRangeForPeptideID == null)
        {
            setRetTimeRangeForPeptideID = con.prepareStatement("UPDATE "+PEPTIDE_TABLE_NAME+
                    " SET startTime = ?, endTime = ? where peptideID = ? ");
        }
        setRetTimeRangeForPeptideID.setFloat(1, (float) startTime);
        setRetTimeRangeForPeptideID.setFloat(2, (float) endTime);
        setRetTimeRangeForPeptideID.setInt(3, id);
        setRetTimeRangeForPeptideID.execute();
    }



    private void setRetTimeRangeForPeptide(String seqKey, double startTime, double endTime) throws SQLException {
        if(setRetTimeRangeForSeqKey == null)
        {
            setRetTimeRangeForSeqKey = con.prepareStatement("UPDATE "+PEPTIDE_TABLE_NAME+
                    " SET startTime = ?, endTime = ? where sequenceCS = ? ");
        }
        setRetTimeRangeForSeqKey.setFloat(1, (float) startTime);
        setRetTimeRangeForSeqKey.setFloat(2, (float) endTime);
        setRetTimeRangeForSeqKey.setString(3, seqKey);
        setRetTimeRangeForSeqKey.execute();
    }


    public void generateRetentionTimeRangesForRedundantSpectraLibrary() throws SQLException {
        List<Integer> peptideIDList = queryPeptideIDS();
        for(int id : peptideIDList)
        {
            boolean isDecoy = isDecoyQuery(id);
            if(!isDecoy)
            {
                List<Double> retTimeList = getAllRetentionTimesForPeptide(id);
                List<Double> filteredList = removeOutlier(retTimeList,0.1);
                double min = Double.MAX_VALUE;
                double max = Double.MIN_VALUE;
                for(double d: filteredList)
                {
                    if(d<min)
                    {
                        min = d;
                    }
                    if(d>max)
                    {
                        max = d;
                    }
                }
                setRetTimeRangeForPeptide(id,min,max);
            }
        }
    }

    public void generateRetTimeRangeForUniqueSpectraLibrary() throws SQLException {
        for(Map.Entry<String,List<Double>> entry: seqKeyRetTimeMap.entrySet())
        {
            String seqKey = entry.getKey();
            List<Double> retTimeList = entry.getValue();
            List<Double> prevList = queryRetentionTimeList(seqKey);
            retTimeList.addAll(prevList);
            List<Double> filteredList = removeOutlier(retTimeList,0.1);
            double min = Double.MAX_VALUE;
            double max = Double.MIN_VALUE;
            for(double d: filteredList)
            {
                if(d<min)
                {
                    min = d;
                }
                if(d>max)
                {
                    max = d;
                }
            }
            setRetTimeRangeForPeptide(seqKey,min,max);
            uploadRetentionTimeList(retTimeList,seqKey);
        }
    }

    public void generateDecoysInMemory() throws SQLException {
        deleteCurrentDecoyDataBase();
        TreeMultimap<Integer,MetaSpectraEntry> map = createMassMap();
        int forwardSize = map.size();
        TreeMultimap<Integer,MetaSpectraEntry> swapMap = generateSwapMap(map);
        uploadDecoyDatabase(swapMap);
        System.out.println(">>decoy data size " +swapMap.size()+"\tforward size: "+forwardSize);
    }

    private void deleteCurrentDecoyDataBase() throws SQLException {
        clearSpectraMetaTable.execute();
    }

    private void uploadDecoyDatabase(TreeMultimap<Integer,MetaSpectraEntry> swapMap) throws SQLException {
        for(Map.Entry<Integer,Collection<MetaSpectraEntry>> entry: swapMap.asMap().entrySet())
        {
            int massKey = entry.getKey();
            Collection<MetaSpectraEntry> specList = entry.getValue();
            for(MetaSpectraEntry currentEntry: specList)
            {
                addDecoyMetaSpectra(currentEntry);
            }
        }
    }


    private TreeMultimap<Integer,MetaSpectraEntry> generateSwapMap(TreeMultimap<Integer,MetaSpectraEntry> forwardMap )
    {
        TreeMultimap<Integer,MetaSpectraEntry> swapMap = TreeMultimap.create();
        NavigableMap<Integer, Collection<MetaSpectraEntry>> navimap = forwardMap.asMap();
        Map.Entry<Integer,Collection<MetaSpectraEntry>> entry;
        while ((entry =navimap.pollFirstEntry())!=null)
        {
            int massKey = entry.getKey();

            /*if(massKey== 2620496)
            {
                System.out.println();
            }*/
            Collection<MetaSpectraEntry> specList = entry.getValue();
            for(MetaSpectraEntry currentEntry: specList)
            {

                int cs = currentEntry.chargeState;
                int diff = cs*DECOY_MASS_SHIFT;
                MetaSpectraEntry swapEntry = findSwapSpectra(massKey, diff,navimap,forwardMap,swapMap);
                float mass_diff = (swapEntry.massKey-currentEntry.massKey)/((float)MZ_KEY_SHIFT);

                MetaSpectraEntry currentDecoy = new MetaSpectraEntry(swapEntry.spectraID,currentEntry.rowId,currentEntry.massKey,currentEntry.chargeState);
                currentDecoy.setDiff(mass_diff/swapEntry.chargeState);

                MetaSpectraEntry swapDecoy = new MetaSpectraEntry(currentEntry.spectraID,swapEntry.rowId, swapEntry.massKey,swapEntry.chargeState );
                swapDecoy.setDiff(-mass_diff/currentEntry.chargeState);
                swapMap.put(currentDecoy.massKey,currentDecoy);
                swapMap.put(swapDecoy.massKey,swapDecoy);
            }
        }
        return swapMap;
    }

    public MetaSpectraEntry  findSwapSpectra(int masskey, int diff, NavigableMap<Integer, Collection<MetaSpectraEntry>> navimap,
                               TreeMultimap<Integer,MetaSpectraEntry> multiMap,TreeMultimap<Integer,MetaSpectraEntry> swapMap )
    {
        int swapPoint = masskey+diff;
       Map.Entry<Integer,Collection<MetaSpectraEntry>> entry = navimap.ceilingEntry(swapPoint);
       if(entry!=null)
       {
           Collection<MetaSpectraEntry> entrylist = entry.getValue();
            if(entrylist.size()>0)
            {
                MetaSpectraEntry result= entrylist.iterator().next();
               // int specID = entrylist.iterator().next().spectraID;
                multiMap.remove(entry.getKey(),result);
                return  result;
            }
       }

        NavigableMap<Integer, Collection<MetaSpectraEntry>> swapNaviMap = swapMap.asMap();
        Map.Entry<Integer,Collection<MetaSpectraEntry>> entryWithDecoy = swapNaviMap.ceilingEntry(swapPoint);

        MetaSpectraEntry result;
        if(entryWithDecoy==null)
        {
            int smallSwappoint = masskey-diff;
            entryWithDecoy = swapNaviMap.floorEntry(smallSwappoint);
            result= entryWithDecoy.getValue().iterator().next();
        }
        else
        {
            result= entryWithDecoy.getValue().iterator().next();
        }
        return result;
    }



    public void fillMetaTable() throws SQLException {
        fillMetaTable.execute();
    }



    public void readDTASelect(String dtaSelectFile) throws JDOMException, SQLException, IOException {
        String path = dtaSelectFile.substring(0,dtaSelectFile.lastIndexOf(File.separatorChar));
        readDTASelect(dtaSelectFile,path);
        if(isUploadAllSpectra())
        {
            uploadAllSpectraToDB(path);
        }
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

    public static Map<String,String> generateDevLineMap(String fasta) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(fasta));
        Map<String,String> devMap = new HashMap<>();
        String line;
        while((line = br.readLine())!=null)
        {
            if(line.startsWith(">"))
            {
                String cleanLine = line.trim().replaceFirst(">","");
                String accession = getSequestLikeAccession(cleanLine);
                devMap.put(accession,cleanLine);
            }

        }
        return devMap;
    }

    public void readDTASelect(String dtaSelectFile,String path) throws IOException, JDOMException, SQLException {
        boolean start = false;
        indexedMap = IndexUtil.createIndexedHt(path,"ms2");
        seqKeyScoreMap = new HashMap<>();
        seqKeyIDMap = new HashMap<>();
        peptideID = queryPeptideCount();
        proteinID = queryProteinCount();
        spectraID = querySpectraCount();
       // initializeLocalMaps();
        int peptideRowLength=0;
        int proteinRowLength =0;
        BufferedReader br = new BufferedReader(new FileReader(dtaSelectFile));
        String line;
        String locus = "";
        String description = "";
        String accession = " ";
        Map<String,String> devLineMap = new HashMap<>();
        List<String> devList = new ArrayList<>();
        List<String> accessionList = new ArrayList<>();

        DTASelectParser parser = new DTASelectParser(dtaSelectFile);
        boolean emptyDevList = false;
        fileNameScanNumSet = new HashSet<>();

        while((line=br.readLine())!=null)
        {
            if(line.endsWith(".fasta"))
            {
                devLineMap = generateDevLineMap(line);
            }
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
                    if(uploadAllSpectra)
                    {
                        String fileNameScanCS = arr[fileNameIndex];
                        fileNameScanNumSet.add(fileNameScanCS);
                    }
                    if(devList.size()>0 && !appendAllSpectra)
                    {
                        processPeptideLine(arr,path,accessionList,devList);
                    }
                    emptyDevList = true;
                }
                else if(arr.length<=proteinRowLength)
                {
                    if(emptyDevList)
                    {
                        emptyDevList = false;
                        accessionList = new ArrayList<>();
                        devList = new ArrayList<>();
                    }
                   locus = arr[locusIndex];
                   description = arr[descriptiveIndex];
                    accession  = getSequestLikeAccession(locus);
                    String devString = devLineMap.get(accession);
                    if(!accession.startsWith("Reverse_") && !accession.startsWith("contaminant"))
                    {
                        accessionList.add(accession);
                        devList.add(description);
                    }
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
        String tempName = filename;
        if(isHeavyFile(filename,path))
        {
            tempName = filename.substring(1,filename.length());
        }
        String spectraKey = path+ tempName + scan;


            IndexedFile ifile = indexedMap.get(tempName);
        List<Float> mzList = new ArrayList<>();
        List<Float> intList = new ArrayList<>();
        //long id = seqKeyIDMap.get(seqKey);
        double retTime = readSpectraFromMS2(ifile,scan,mzList,intList);
        List<Double> retTimeList = seqKeyRetTimeMap.getOrDefault(seqKey,new ArrayList<>());
        retTimeList.add(retTime);
        seqKeyRetTimeMap.put(seqKey,retTimeList);

        updateEntry(seqKey,mass,mzList.size(),0,(float)retTime,0,0,filename,
                score,0,deltacn,scan);

       // if(!fileNameScanNumSet.contains(spectraKey))
        {
            updateSpectra(id,mzList,intList,mass,(float)retTime, filename,scan);
          //  fileNameScanNumSet.add(spectraKey);
        }

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
        float retTime = (float)readSpectraFromMS2(ifile,scan,mzList,intList);


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

    public void uploadAllSpectraToDB(String path) throws IOException, SQLException {
        File dir = new File(path);
        File [] ms2Arr = dir.listFiles(file -> file.getName().endsWith(".ms2"));
        Set<String> fileSet = new HashSet<>();
        for(File ms2 : ms2Arr)
        {
            String fileName = ms2.getName();
            if(fileName.startsWith("H")
                    && new File(path+File.separatorChar + fileName.substring(1)).exists())
            {
                continue;
            }
            SpectrumReader reader = new SpectrumReader(ms2.getAbsolutePath(), "ms2");

            for(Iterator<PeakList> itr = reader.getSpectra(); itr.hasNext(); )
            {
                PeakList peakList = itr.next();
                float retTime =(float) peakList.getRetentionTime();

                int scan = peakList.getLoscan();
                List<Zline> zlineToUpload = new ArrayList<>();

                String fileNameScan = fileName+"."+scan+"."+scan+".";
                for(Iterator<Zline> ztr = peakList.getZlines(); ztr.hasNext();)
                {
                    Zline zline = ztr.next();
                    int cs  = zline.getChargeState();
                    String fileNameScanCS = fileNameScan + cs;
                    if(!fileNameScanNumSet.contains(fileNameScanCS))
                    {
                        zlineToUpload.add(zline);
                    }

                }

                if(zlineToUpload.size()>0)
                {
                    List<Float> massList = new ArrayList<>();
                    List<Float> intList =  new ArrayList<>();
                    for(Peak peak: peakList.getPeaksList())
                    {
                        massList.add((float) peak.getM2z());
                        intList.add((float)peak.getIntensity());
                    }
                    for(Zline zline: zlineToUpload)
                    {
                        float prcMass = (float) zline.getM2z();
                        int cs = zline.getChargeState();
                        uploadSpectraWithoutPeptide(scan,cs,prcMass,retTime,fileName, massList, intList);
                    }
                }

            }

            reader.closeDataFile();
        }
    }

    public static String getWeightCode(String fileName, String filePath)
    {
        String weightCode = "";
        weightCode = isHeavyFile(fileName,filePath) ? "H" : weightCode;
        weightCode = isMediumFile(fileName,filePath) ? "M" : weightCode;
        return  weightCode;
    }


    private void processPeptideLine(String[] peptideArr,String path, List<String> accessionList,
                                    List<String> descriptionList) throws IOException, SQLException {
        float xcorr = Float.parseFloat(peptideArr[xcorrIndex]);
        float deltaCN = Float.parseFloat(peptideArr[deltaCNIndex]);
        String sequence = peptideArr[seqIndex];
        float theoreticalMass = Float.parseFloat(peptideArr[theoreticalMassIndex]);
        float experimentalMass = Float.parseFloat(peptideArr[theoreticalMassIndex-1]);
        float mass = theoreticalMass;
        String fileString = peptideArr[fileNameIndex];

        String regex = "\\.";
        String weightCode = "";
        String [] fileArr = fileString.split(regex);

        String fileName=  fileArr[0]+".ms2";
        weightCode = isHeavyFile(fileName,path) ? "H" : weightCode;
        weightCode = isMediumFile(fileName,path) ? "M" : weightCode;
        int cs = Integer.parseInt(fileArr[3]);
        int scan = Integer.parseInt(fileArr[2]);
        long peptideID=-1;
        String seqKey = removePrePostFix(sequence)+cs +weightCode;


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
            if(isHeavyFile(fpath,path))
            {
                fpath = fileName.substring(1,fileName.length());
            }
            if(new File(path+File.separator+fpath).exists())
            {
                IndexedFile ifile = indexedMap.get(fpath);

                List<Float> mzList = new ArrayList<>();
                List<Float> intList = new ArrayList<>();
                double retTime = readSpectraFromMS2(ifile,scan,mzList,intList);
                List<Double> retTimeList = seqKeyRetTimeMap.getOrDefault(seqKey,new ArrayList<>());
                retTimeList.add(retTime);
                seqKeyRetTimeMap.put(seqKey,retTimeList);

                peptideID = insertEntry(sequence,mass,cs,mzList.size(),(float) retTime,0,0,fileName,
                        xcorr,0,deltaCN,scan,seqKey);
                seqKeyIDMap.put(seqKey,peptideID);
                proteinSet = new TLongHashSet();
                seqKeyProteinSetMap.put(seqKey, proteinSet);
                seqKeyScoreMap.put(seqKey,xcorr);

                insertSpectra(peptideID,mass,mzList,intList, (float)retTime, fileName, scan);

            }
            else
            {
                return;
            }

        }
        /*
          check if new score is higher than library score
                upload if true
        */
        else if(xcorr>libraryScore && uploadBestScoringXcorr) {
            seqKeyScoreMap.put(seqKey,xcorr);
            if(new File(path+File.separator+fileName).exists())
            {
                readAndUpdateSpectra(seqKey, fileName, path, xcorr, mass, deltaCN, scan, peptideID);
            }
            else
            {
                return;
            }
        }
        else if(!uploadBestScoringXcorr)
        {
            String key = fileString.trim();

            String fpath = fileName;
            if(isHeavyFile(fpath,path))
            {
                fpath = fileName.substring(1,fileName.length());
            }

            if(!fileStringSet.contains(key) && new File(path+File.separator+fpath).exists())
            {
                fileStringSet.add(key);

                IndexedFile ifile = indexedMap.get(fpath);
                List<Float> mzList = new ArrayList<>();
                List<Float> intList = new ArrayList<>();

                String spectraKey = path+ fileName + scan;
                if(!fileNameScanNumSet.contains(spectraKey))
                {
                    float retTime = (float)readSpectraFromMS2(ifile,scan,mzList,intList);
                    insertSpectra(peptideID,mass,mzList,intList,retTime, fileName, scan);
                }
            }

        }
        /*
               Check if accession exists in local memory
                    if not check if exists in library
                        if not add protein to protein table
                            flag to add to peptide
                    else
                        add to local memory
         */
        for(int i=0; i<accessionList.size(); i++)
        {
            String accession = accessionList.get(i);
            String description = descriptionList.get(i);
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
                    proteinID = insertProtein(accession,description);
                }
            }
            else proteinID = locusProteinIDMap.get(seqKey);
            if(addProteinToPeptide || !proteinSet.contains(proteinID))
            {
                proteinSet.add(proteinID);
                //System.out.println(peptideID+ " " +proteinID);
                insertProteinToPeptideIndex(peptideID,proteinID);
            }

        }

    }

    private double readSpectraFromMS2(IndexedFile ifile, int scan,List<Float> mzList, List<Float> intList)
            throws IOException {
        if(ifile== null )
            System.out.println("WHOOPS!");
        //System.out.println(">>>>"+ifile.getFileName()+"\t"+scan);
        long pos = ifile.getPositionByScan(scan);
        RandomAccessFile rfile = new RandomAccessFile(ifile.getFileName(),"r");
        rfile.seek(pos);
        double retTime = -1;
        String line=rfile.readLine();

        while((line=rfile.readLine())!=null && !line.startsWith("S")){
            if(line.startsWith("I\tRetTime\t"))
            {
                retTime = Double.parseDouble( line.split("\t")[2]);

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
       queryEntry(mzStart,mzEnd,result);
       if(result.size()>0)
       {
           querySpectra(mzStart,mzEnd,result);
       }
        return result;
    }





    public static void main(String [] args) throws SQLException, IOException, JDOMException {
        //testLibCreation();
        //testReadingDTASelect();
        //testReadLibrary();
//        testReadLib2();
        String libAddress = args[0];
        System.out.println(">>library is "+libAddress);
        File f = new File(libAddress);

        boolean justGenerateDecoys =false;
        boolean uniqueSpectra = true;
        boolean keepOldDatabase  = false;
        boolean uploadAllSpectra = false;
        boolean appendAllSpectra = false;
        boolean skipDecoys = false;
        boolean skipDTA =false;
        for(String s: args  )
        {
            String slower = s.toLowerCase();
            if( slower.equals("--decoys"))
            {
                skipDTA = true;
                justGenerateDecoys = true;
                keepOldDatabase = true;
            }
            else if(slower.equals("--redundant-spectra"))
            {
                uniqueSpectra = false;
            }
            else if(slower.equals("--append-old-db"))
            {
                keepOldDatabase = true;
            }
            else if(slower.equals("--upload-all-spectra"))
            {
                uploadAllSpectra = true;
            }
            else if(slower.endsWith("--append-all-spectra"))
            {
                appendAllSpectra = true;
                uploadAllSpectra = true;
                keepOldDatabase = true;
            }
            else if(slower.equals("--skip-decoys"))
            {
                skipDecoys = true;
            }
            else if(slower.equals("--skip-dta"))
            {
                skipDTA = true;
                keepOldDatabase = true;
            }
        }
        if (!keepOldDatabase && f.exists())
            f.delete();

        LibraryIndexer libraryIndexer = new LibraryIndexer(libAddress);
        libraryIndexer.setUploadAllSpectra(uploadAllSpectra);
        libraryIndexer.setUploadBestScoringXcorr(uniqueSpectra);
        libraryIndexer.setAppendAllSpectra(appendAllSpectra);
        if(!skipDTA)
        {
            for(int i=1; i<args.length; i++)
            {
                if(!args[i].startsWith("--"))
                {
                    System.out.println("reading "+args[i]);
                    libraryIndexer.readDTASelect(args[i]);

                }

            }
        }
        if(uniqueSpectra)
        {
           libraryIndexer.generateRetTimeRangeForUniqueSpectraLibrary();
        }
        else
        {
            libraryIndexer.generateRetentionTimeRangesForRedundantSpectraLibrary();
        }

        if(!skipDecoys)
        {
            libraryIndexer.generateDecoysInMemory();
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

    public int getNextHighestMassIndex(int massIndex) throws SQLException {
        getNextHighestMassIndex.setInt(1,massIndex);
        ResultSet rs = getNextHighestMassIndex.executeQuery();
        int nextIndex = 0;
        while( rs.next() )
        {
            nextIndex = rs.getInt(1);
        }
        return nextIndex;
    }

    public void createSpectraMetaTable() throws SQLException {
        createSpectraMetaTable.execute();
    }

    public int queryDecoyCount() throws SQLException {
        ResultSet set = getNumberOfDecoysStatement.executeQuery();
        if(set.next())
        {
            return set.getInt(1);
        }
        return -1;
    }



    public List<LibrarySpectra> getNPeptides(int massIndex, int n) throws SQLException {

        getNPeptidesStatement.setInt(1,massIndex);
        getNPeptidesStatement.setInt(2,n);
        final ResultSet rs = getNPeptidesStatement.executeQuery();
        int i=0;
        List<LibrarySpectra> libList = new ArrayList<>();

        long prevId =-1;
        LibrarySpectra prevsSpectra = null;
        while (rs.next()) {
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
                int massKey = rs.getInt(3);
                float xcorr = rs.getFloat(12);
                float mass = rs.getFloat(4);
                int cs = rs.getInt(5);
                float retTime = rs.getFloat(8);
                float deltaCN = rs.getFloat(14);
                String key = rs.getString(15);
                String fileName = rs.getString(11);
                int scan = rs.getInt(16);

                LibrarySpectra librarySpectra = new LibrarySpectra(seq,cs,mass,retTime, xcorr, deltaCN,key,fileName,scan,id, accession, proteinDescription, massKey);
                libList.add(librarySpectra);
                prevsSpectra =librarySpectra;
            }
            else prevsSpectra.addProtein(accession,proteinDescription);
            prevId=id;
            i++;
        }
        return  libList;
    }

    public List<MetaSpectraEntry> getMetaSpectra(int massIndex, int num) throws SQLException {

        List<MetaSpectraEntry> entryList = new ArrayList<>();
        getNMetaAllSpectraStatement.setInt(1, massIndex);
        //getNMetaForwardSpectraStatement.setInt(2, mod);
        getNMetaAllSpectraStatement.setInt(2, num);
        final ResultSet rs = getNMetaAllSpectraStatement.executeQuery();

        while (rs.next()) {
            int rowId = rs.getInt(1);
            int spectraID = rs.getInt(2);
            int massKey = rs.getInt(3);
            int chargeState = rs.getInt(4);
            MetaSpectraEntry entry = new MetaSpectraEntry(rowId,spectraID,massKey, chargeState);
            entryList.add(entry);

        }
        rs.close();
        return entryList;
    }

    public List<MetaSpectraEntry> getForwardMetaSpectra(int massIndex, int num) throws SQLException {

        List<MetaSpectraEntry> entryList = new ArrayList<>();
        getNMetaForwardSpectraStatement.setInt(1, massIndex);
        //getNMetaForwardSpectraStatement.setInt(2, mod);
        getNMetaForwardSpectraStatement.setInt(2, num);
        final ResultSet rs = getNMetaForwardSpectraStatement.executeQuery();

        while (rs.next()) {
            int rowId = rs.getInt(1);
            int spectraID = rs.getInt(2);
            int massKey = rs.getInt(3);
            int chargeState = rs.getInt(4);
            MetaSpectraEntry entry = new MetaSpectraEntry(rowId,spectraID,massKey, chargeState);
            entryList.add(entry);

        }
        rs.close();
        return entryList;
    }

    public void addDecoyMetaSpectra(MetaSpectraEntry decoyEntry) throws SQLException {
        addDecoyMetaSpectra(decoyEntry.spectraID, decoyEntry.massKey, decoyEntry.chargeState, decoyEntry.diff );
    }


    public void addDecoyMetaSpectra(int spectraID, int massKey, int chargeState, float diff) throws SQLException {
        insertSpectraMetaTable.setInt(1,spectraID);
        insertSpectraMetaTable.setInt(2,massKey);
        //   insertSpectraMetaTable.setInt(3,decoyEntry.rowId);
        insertSpectraMetaTable.setInt(3,chargeState);
        insertSpectraMetaTable.setInt(4,1);
        insertSpectraMetaTable.setInt(5,0);
        insertSpectraMetaTable.setFloat(6,diff);
        insertSpectraMetaTable.execute();
    }



    public static class MetaSpectraEntry implements Comparable<MetaSpectraEntry>
    {
        public final int spectraID;
        public final int rowId;
        public final int massKey;
        public final int chargeState;
        private int decoyId = -1;
        private float diff = 0;

        public MetaSpectraEntry(int spectraID, int rowId, int massKey, int chargeState) {
            this.spectraID = spectraID;
            this.rowId = rowId;
            this.massKey = massKey;
            this.chargeState = chargeState;
        }

        @Override
        public int compareTo(MetaSpectraEntry metaSpectraEntry) {
            return rowId-metaSpectraEntry.rowId;
        }

        public float getDiff() {
            return diff;
        }

        public void setDiff(float diff) {
            this.diff = diff;
        }
    }

    private void updateHasDecoyState(int massStart, int massEnd) throws SQLException {
        updateHasDecoyStateStatement.setInt(1,massStart);
        updateHasDecoyStateStatement.setInt(2,massEnd);
        //updateHasDecoyStateStatement.setInt(3,mod);
        updateHasDecoyStateStatement.execute();
    }

    public TreeMultimap<Integer, MetaSpectraEntry> createMassMap() throws SQLException {
        final ResultSet rs = getAllSpectraInformation.executeQuery();
        TreeMultimap<Integer,MetaSpectraEntry> treeMultimap = TreeMultimap.create();
        int i=0;
        while(rs.next())
        {
            int rowId = rs.getInt(1);
            int massKey = rs.getInt(2);
         /*   if(massKey==2620496)
            {
                System.out.println();
            }*/
            int chargeState = rs.getInt(3);
            MetaSpectraEntry entry = new MetaSpectraEntry(rowId,i++,massKey,chargeState);
            treeMultimap.put(massKey,entry);
        }
        rs.close();
        return treeMultimap;
    }

    public boolean isUploadBestScoringXcorr() {
        return uploadBestScoringXcorr;
    }

    public void setUploadBestScoringXcorr(boolean uploadBestScoringXcorr) {
        this.uploadBestScoringXcorr = uploadBestScoringXcorr;
    }

    public static boolean isHeavyFile(String file, String path) {

        if(!file.startsWith("H"))
            return false;


        String lightfile = path + "/" + file.substring(1, file.length());

        File lf = new File(lightfile);

        if(lf.exists())
            return true;
        else
            return false;

    }


    public static boolean isMediumFile(String file, String path) {

        if(!file.startsWith("M"))
            return false;


        String lightfile = path + "/" + file.substring(1, file.length());

        File lf = new File(lightfile);

        if(lf.exists())
            return true;
        else
            return false;

    }

    private void uploadRetentionTimeList(List<Double> retentionTimeList, String seqKey) throws SQLException {
        StringBuilder sb = new StringBuilder();
        for(double d: retentionTimeList)
        {
            sb.append(d).append(";");
        }
        if(setRetentionTimeList==null)
        {
            setRetentionTimeList = con.prepareStatement("UPDATE "+PEPTIDE_TABLE_NAME
                    +" SET retentionTimeList = ? WHERE sequenceCS = ?");
        }

        setRetentionTimeList.setString(1,sb.toString());
        setRetentionTimeList.setString(2,seqKey);
        setRetentionTimeList.execute();
    }

    private List<Double> queryRetentionTimeList(String seqKey) throws SQLException {
        if(getRetentionTimeList == null)
        {
            getRetentionTimeList = con.prepareStatement("SELECT retentionTimeList FROM "+PEPTIDE_TABLE_NAME+
                    " WHERE sequenceCS = ?");
        }
        getRetentionTimeList.setString(1,seqKey);
        ResultSet rs = getRetentionTimeList.executeQuery();
        String listStr = rs.getString(1);
        List<Double> retentionTimeList = new ArrayList<>();

        if(listStr!=null && listStr.length()>0)
        {
            String [] arr  = listStr.split(";");
            for(String s: arr)
            {
                double d = Double.parseDouble(s);
                retentionTimeList.add(d);
            }
        }

        return  retentionTimeList;
    }

    public void setUploadAllSpectra(boolean uploadAllSpectra) {
        this.uploadAllSpectra = uploadAllSpectra;
    }

    public boolean isUploadAllSpectra() {
        return uploadAllSpectra;
    }

    public void setAppendAllSpectra(boolean appendAllSpectra) {
        this.appendAllSpectra = appendAllSpectra;
    }

    public void uploadSpectraWithoutPeptide(int scanNumber, int cs, float prcMass, float retTime, String fileName,
                                            List<Float> mzList, List<Float> intList) throws IOException, SQLException {
        String code = Integer.toString((int)peptideID);

        String sequence = "N.OTAPEPTID.E"+(code);
        String seqKey = sequence+cs;
        long oldId = insertEntry("",prcMass,cs,mzList.size(),(float) retTime,0,0,fileName,
                0,0,0,scanNumber,seqKey, true);
        insertSpectra(oldId, prcMass, mzList, intList, retTime, fileName, scanNumber, true);
        int massKey = (int)(prcMass*MZ_KEY_SHIFT);
       // addDecoyMetaSpectra((int)spectraID, massKey, cs, 0  );

        //peptideID++;

    }

    private PreparedStatement isDecoyStatement = null;
    public boolean isDecoyQuery(int id) throws SQLException {
        if(isDecoyStatement==null)
        {
            isDecoyStatement = con.prepareStatement("select isDecoy from peptidetable where peptideid = ?");
        }
        isDecoyStatement.setInt(1,id);
        ResultSet rs = isDecoyStatement.executeQuery();
        boolean result =false;
        while(rs.next())
        {
            result = rs.getInt(1)==1 ;
        }
        return result;
    }



}
