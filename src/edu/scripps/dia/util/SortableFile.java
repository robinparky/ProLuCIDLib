
/*
* Copyright (c) 2008 Proteomics Integrated Solutions. All rights reserved.  
*/


/**
 *
 * @author Tao Xu
 * @email taoxu@integratedproteomics.com
 * Created on Nov 11, 2010
 * $Revision
 * $Date
 */

package edu.scripps.dia.util;

import java.io.File;

public class SortableFile implements Comparable {

    private String fileName;
    private long fileSize;
    public SortableFile(String file) {
        fileName = file;
        fileSize = new File(file).length();

    }

    public String getFileName() {
        return fileName;
    } 
    public long getFileSize() {
        return fileSize;
    }
 
    public int compareTo(Object obj) {
        SortableFile sf = (SortableFile) obj;

        long diff = fileSize - sf.fileSize;
        if( diff == 0 ) {
            return 0;
        } else if(diff < 0) {
            return -1;

        } else {
            return 1;
        }
    }
} 
