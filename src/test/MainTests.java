package test;

import org.testng.TestNG;

/**
 * Created by Titus Jung titusj@scripps.edu on 4/9/18.
 */
public class MainTests {


    public static void main(String[] args) {

        TestNG testNG = new TestNG();
        Class[] classes = new Class[]{test.LibraryIndexerTest.class};
        testNG.setTestClasses(classes);
        testNG.run();
    }


}
