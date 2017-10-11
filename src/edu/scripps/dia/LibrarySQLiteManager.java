package edu.scripps.dia;

import gnu.trove.TFloatArrayList;
import gnu.trove.TIntArrayList;
import org.sqlite.SQLiteConfig;

import java.sql.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by yateslab on 8/7/17.
 */
public class LibrarySQLiteManager {

    protected Connection con;
    private String dbPath;
    private PreparedStatement addSeqStatement;
    private PreparedStatement retrieveSeqStatement;
    private List<String> result = new ArrayList<>();


    public static void main(String [] args) throws SQLException {
        String dbPath = args[0];
        LibrarySQLiteManager manager = new LibrarySQLiteManager(dbPath);
        manager.retrieveSpectra(970_000,980_000);
    }


    public LibrarySQLiteManager(String dbPath) throws SQLException {
        this.dbPath = dbPath;
        setUp();
    }

    private void setUp() throws SQLException {

        SQLiteConfig config = new SQLiteConfig();
        //optimize for multiple connections that can share data structures
        config.setSharedCache(true);
        //config.setCacheSize(cacheSize);
        //config.setPageSize(pageSize);
        config.setJournalMode(SQLiteConfig.JournalMode.OFF);
        config.enableEmptyResultCallBacks(false);
        config.enableCountChanges(false);
        config.enableFullSync(false);
        config.enableRecursiveTriggers(false);
        config.setLockingMode(SQLiteConfig.LockingMode.NORMAL);
        config.setSynchronous(SQLiteConfig.SynchronousMode.OFF); //TODO may be dangerous on some systems to have off


        con = DriverManager.getConnection("jdbc:sqlite:" + dbPath,config.toProperties());

        con.prepareStatement("CREATE TABLE IF NOT EXISTS library_sequence "
                + "(precursor_mass_key INTEGER , "+"charge_state INTEGER,"
                + "sequence VARCHAR,"+"massList VARCHAR,"+"intensityList VARCHAR"
                + ");").execute();
        con.prepareStatement("CREATE INDEX IF NOT EXISTS precursor_mass_key_index_dsc ON "
                + getIndexTableName() +" (precursor_mass_key DESC);").execute();

        addSeqStatement = con.prepareStatement(
                "INSERT INTO " + getIndexTableName() + " (precursor_mass_key, charge_state,sequence, massList, intensityList) "
                        + "VALUES (?, ?, ?, ?, ?);");
        retrieveSeqStatement = con.prepareStatement("SELECT precursor_mass_key, sequence, massList, intensityList "
                + "FROM " + getIndexTableName() + " "
                + "WHERE precursor_mass_key BETWEEN ? AND ?;");
    }

    public void loadSpectra(String seq, int mass, int cs, TFloatArrayList massList, TIntArrayList intensityList) throws SQLException {
        addSeqStatement.setInt(1,mass);
        addSeqStatement.setInt(2,cs);
        addSeqStatement.setString(3,seq);
        addSeqStatement.setString(4,massList.toString());
        addSeqStatement.setString(5,intensityList.toString());

        addSeqStatement.executeUpdate();
    }

    public void retrieveSpectra(int startRange, int endRange) throws SQLException {
        retrieveSeqStatement.setInt(1,startRange);
        retrieveSeqStatement.setInt(2,endRange);
        ResultSet rs = retrieveSeqStatement.executeQuery();
        long start = System.currentTimeMillis();
        while(rs.next())
        {
            int mass = rs.getInt(1);
            String seq =rs.getString(2);
            String massStringList=rs.getString(3);
            String intensityStringList=rs.getString(4);
            String [] massStrArr = massStringList.substring(1,massStringList.length()-2).split(",");
            System.out.println("mass "+mass);
           // System.out.println("massString "+ Arrays.toString(massStrArr));
            //System.out.println("intensityStringList "+intensityStringList);
        }
        long duration = System.currentTimeMillis() - start;
        System.out.println(duration);
    }

    private String retrieveSequence(ResultSet rs) throws SQLException {
        return rs.getString(2);
    }

    private TFloatArrayList retrieveMassList(ResultSet rs) throws SQLException {
        String massStringList=rs.getString(3);
        String [] massStrArr = massStringList.substring(1,massStringList.length()-2).split(",");
        TFloatArrayList result = new TFloatArrayList();
        for(String s:massStrArr)
        {
            result.add(Float.parseFloat(s));
        }
        return result;
    }

    private TIntArrayList retrieveIntensityList(ResultSet rs) throws SQLException {
        String intensityStr=rs.getString(4);
        String [] intensityStrArray = intensityStr.substring(1,intensityStr.length()-2).split(",");
        TIntArrayList result = new TIntArrayList();
        for(String s:intensityStrArray)
        {
            result.add(Integer.parseInt(s));
        }
        return result;
    }

    private int[] retrieveMassArray(ResultSet rs, int massArraySize) throws SQLException {
        String intensityStr=rs.getString(4);
        String [] intensityStrArray = intensityStr.substring(1,intensityStr.length()-2).split(",");
        String massStringList=rs.getString(3);
        String [] massStrArr = massStringList.substring(1,massStringList.length()-2).split(",");
        int [] massArray = new int[massArraySize];

        for(int i=0; i<intensityStrArray.length; i++)
        {
            int mass = (int) (Float.parseFloat( massStrArr[i])*1000);
            int intensity = Integer.parseInt(intensityStrArray[i]);
            if(mass>=massArraySize) break;
            massArray[mass] = intensity;
        }
        return massArray;

    }

    private int retrieveMass(ResultSet rs) throws SQLException {
        return rs.getInt(1);
    }

    private QueryResult createQueryResult(ResultSet rs, int massArraySize) throws SQLException {
        String seq = retrieveSequence(rs);
        float mass = retrieveMass(rs);
        int[] massArray = retrieveMassArray(rs,massArraySize);
        return new QueryResult(seq,mass,massArray);
    }

    public List<QueryResult> retrieveSpectra(int startRange, int endRange, int massArraySize) throws SQLException {
        retrieveSeqStatement.setInt(1,startRange);
        retrieveSeqStatement.setInt(2,endRange);
        ResultSet rs = retrieveSeqStatement.executeQuery();
        List<QueryResult> queryResultList = new ArrayList<>();
        while(rs.next())
        {
            queryResultList.add(createQueryResult(rs,massArraySize));
        }
        return queryResultList;
    }




    public String getIndexTableName() {
        return "library_sequence";
    }
}
