package edu.scripps.dia.models;

/**
 * Created by yateslab on 7/19/17.
 */
public class FileScan {

    public String file;
    public String scan;
    public FileScan(String file, String Scan)
    {
        this.file = file;
        this.scan = Scan;
    }

    @Override
    public int hashCode() {
        return file.hashCode()+scan.hashCode();
    }
}
