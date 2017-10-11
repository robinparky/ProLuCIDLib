/**
 * @file FileFormatUnknownException.java
 * This is the source file for edu.scripps.pms.util.spectrum.FileFormatUnknowException
 * @author Tao Xu
 * @date $Date: 2007/02/26 23:47:37 $
 */



package edu.scripps.dia.util;

import java.io.IOException;

public class FileFormatUnknownException extends IOException {
    
    private static String messageStr = "Unknow file format";

    public FileFormatUnknownException() {

        super(messageStr);
    }
}
