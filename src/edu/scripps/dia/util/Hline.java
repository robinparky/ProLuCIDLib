package edu.scripps.dia.util;

import java.util.Iterator;
import java.util.List;

public class Hline {
    private String extractor;
    private String extractorVersion;
    private String comments;
    private String extractorOptions;
    private String creationDate;


    public Hline(List hlines) {
        init(hlines.iterator());
    }

    public Hline(Iterator<String> hlines) {
        init(hlines);
    }
    public void init(Iterator<String> itr)
    {
        String eachLine;
        String[] strArr;

        while (itr.hasNext())
        {
            eachLine = itr.next().toString();
            strArr = eachLine.split("\t");

            if(null == strArr[1])
                continue;

            if(strArr[1].endsWith("tor"))
                setExtractor(strArr[2]);
            else if(strArr[1].endsWith("sion"))
                setExtractorVersion(strArr[2]);
            else if(strArr[1].endsWith("ents"))
                setComments(strArr[2]);
            else if(strArr[1].endsWith("tions"))
                setExtractorOptions(strArr[2]);
            else if(strArr[1].endsWith("Date"))
                setCreationDate(strArr[2]);
        }
    }

    public void setExtractorOptions(String extractorOptions) {
        this.extractorOptions = extractorOptions;
    }

    public void setComments(String comments) {
        this.comments = comments;
    }

    public void setExtractorVersion(String extractorVersion) {
        this.extractorVersion = extractorVersion;
    }

    public void setExtractor(String extractor) {
        this.extractor = extractor;
    }

    public void setCreationDate(String creationDate) {
        this.creationDate = creationDate;
    }

    public String getExtractor() {
        return extractor;
    }

    public String getExtractorVersion() {
        return extractorVersion;
    }

    public String getComments() {
        return comments;
    }

    public String getExtractorOptions() {
        return extractorOptions;
    }

    public String getCreationDate() {
        return creationDate;
    }
}
