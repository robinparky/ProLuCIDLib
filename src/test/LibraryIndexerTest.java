package test;

import edu.scripps.dia.LibraryIndexer;
import edu.scripps.dia.LibrarySpectra;
import gnu.trove.TFloatArrayList;
import org.jdom.JDOMException;
import org.testng.AssertJUnit;
import org.testng.ITestResult;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;


/**
 * Created by Titus Jung titusj@scripps.edu on 4/9/18.
 */
public class LibraryIndexerTest {
    public static final String ms2Address = "testdata/0604TestDTASelect/small.ms2";

    public static final String dbAddress = "testdata/0604TestDTASelect/small.lib";
    public static final String dtaAddress = "testdata/0604TestDTASelect/small.txt";
    public static final String dta2Address = "testdata/0604TestDTASelect/small2.txt";
    public static final String address = "testdata/0604TestDTASelect/";
    public static Set<String> peptideKeysFromDTASet = new HashSet<>();
    public static List<String> peptideList = new ArrayList<>();
    public static LibraryIndexer libraryIndexer;


    @BeforeClass
    public static void StartClean() throws SQLException, IOException, JDOMException {
        deleteLibrary();
        try {
            libraryIndexer = new LibraryIndexer(dbAddress);
            libraryIndexer.readDTASelect(dtaAddress, address);
            //libraryIndexer.readDTASelect(dta2Address, address);

        }
        catch(Exception e){
            e.printStackTrace();
        }
    }


    @Test
    public static void createIndexBasicTest() {
        //AssertJUnit.assertTrue(!new File(dbAddress).exists());
      //  LibraryIndexer libraryIndexer = new LibraryIndexer(dbAddress);
       // AssertJUnit.assertTrue(new File(dbAddress).exists());
        AssertJUnit.assertTrue(new File(dbAddress).exists());
    }

    @Test
    public static void checkIfValidPeptideID() throws SQLException {
        List<Integer> peptideIDTable = new ArrayList<>();
        List<Integer> proteinIDTable = new ArrayList<>();
        libraryIndexer.queryPeptideProteinTable(peptideIDTable,proteinIDTable);
       // Collections.sort(peptideIDTable);
        for(int i: peptideIDTable)
        {
            AssertJUnit.assertTrue(i>=0);
        }
    }

    @Test
    public static void checkIfPeptideProteinDuplicates() throws SQLException {
        List<Integer> peptideIDTable = new ArrayList<>();
        List<Integer> proteinIDTable = new ArrayList<>();
        libraryIndexer.queryPeptideProteinTable(peptideIDTable,proteinIDTable);
        Map<Integer,Set<Integer>> peptideMap = new HashMap<>();
        for(int i=0; i<peptideIDTable.size(); i++)
        {
            int pepID = peptideIDTable.get(i);
            int protID = proteinIDTable.get(i);
            Set<Integer> set = peptideMap.get(pepID);
            if(set!=null)
            {
                AssertJUnit.assertTrue(!set.contains(protID));
                set.add(protID);
            }
            else
            {
                set = new HashSet<>();
                set.add(protID);
                peptideMap.put(pepID,set);
            }
        }
    }
    @Test
    public static void checkPeptideIDsValid() throws SQLException {
        List<Integer> ids = libraryIndexer.queryPeptideIDS();
        for(int i: ids)
        {
            AssertJUnit.assertTrue(i>=0);
        }
    }

    @Test
    public static void checkPeptideIDsUnique() throws SQLException {
        List<Integer> ids = libraryIndexer.queryPeptideIDS();
        Set<Integer> set = new HashSet<>();
        for(int i: ids)
        {
            AssertJUnit.assertTrue(!set.contains(i));
            set.add(i);
        }
    }
    @Test
    public static void checkPeptideKeysUnique() throws SQLException {
        List<String> peptideKeysList = libraryIndexer.queryPeptideKeys();
        Set<String> set = new HashSet<>();
        for(String s: peptideKeysList)
        {
            AssertJUnit.assertTrue(!set.contains(s));
            set.add(s);
        }
    }
    @Test
    public static void checkPeptideEntries() throws SQLException, IOException {
        float start =1800;
        float end = 1850;
        List<LibrarySpectra> libList =libraryIndexer.querySpectra(start,end);
        boolean containsSequence1 = false;
        boolean containsSequence2 = false;
        for(LibrarySpectra spectra: libList)
        {
            if(spectra.sequence.equals("K.GLVAVITGGASGLGLATAER.L")
                    && spectra.chargeState==2
                    && Float.compare(spectra.score,5.6305f)==0
                    && spectra.scan==54202
                    )
            {
                TFloatArrayList mzList = spectra.getMzList();
                TFloatArrayList intList = spectra.getIntensityList();
                //spectra.getAccessionList().contains("IPI:IPI00218343.4|SWISS-PROT:Q9BQE3|TREM");
                if(Float.compare(mzList.get(0),101.071f)==0 && Float.compare(mzList.get(mzList.size()-1),1472.786f)==0
                        && Float.compare(intList.get(0),7626.5f)==0 && Float.compare(intList.get(intList.size()-1),3616.4f)==0 &&
                        spectra.getAccessionList().contains("IPI:IPI00218343.4|SWISS-PROT:Q9BQE3|TREM") &&
                        spectra.getAccessionList().contains("sp|Q99714|HCD2_HUMAN"))
                containsSequence1 = true;
            }
            if(spectra.sequence.equals("K.GLVAVITGGASGLGLATAER.L")
                    && spectra.chargeState==3
                    && Float.compare(spectra.score,4.0581f)==0
                    && spectra.scan==54212
                    )
            {
                TFloatArrayList mzList = spectra.getMzList();
                TFloatArrayList intList = spectra.getIntensityList();

                if(Float.compare(mzList.get(0),101.0712f)==0 && Float.compare(mzList.get(mzList.size()-1),1262.6685f)==0
                        && Float.compare(intList.get(0),11955.1f)==0 && Float.compare(intList.get(intList.size()-1),2043.3f)==0)
                containsSequence2 = true;
            }
        }
        AssertJUnit.assertTrue(containsSequence1&&containsSequence2);

    }


    @AfterMethod
    public void tearDown(ITestResult result) {
        if (result.getStatus() == ITestResult.FAILURE) {
            System.out.println("failed at " +result.getName());
            //your screenshooting code goes here
        }
    }

    @AfterClass
    public void clean()
    {
        deleteLibrary();
    }

    public static void deleteLibrary()
    {
        File db = new File(dbAddress);
        if(db.exists()) db.delete();
    }

}
