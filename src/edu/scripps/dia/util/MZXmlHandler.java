package edu.scripps.dia.util;

import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import org.systemsbiology.jrap.Base64;

import java.nio.ByteBuffer;

/**
 * @author Robin Park
 * @version $Id: MZXmlHandler.java,v 1.3 2009/10/26 22:58:23 taoxu Exp $
 */
public class MZXmlHandler 
{

    public static byte [] floatTobyte(float num)
    {
	ByteBuffer buf = ByteBuffer.allocate(4);
	buf.putFloat(num);
	return buf.array();
    }



    public static PeakList decode64(String input)
    {
	byte[] byteArr = Base64.decode(input);

	int floatBytes = 64 / 8;
	int fieldIndex = 2;
	int i;
	double mass=-1;
	double intensity=-1;
	PeakList pList = new PeakList();

	if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	for (i = 0; i < byteArr.length - floatBytes + 8; i += floatBytes)
	{
	    long intBits = 0;
	    intBits |= (((long) byteArr[i]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 1]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 2]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 3]) & 0xff);
	    // Must be in IEEE 754 encoding!
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 4]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 5]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 6]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 7]) & 0xff);

	    if (fieldIndex == 2)
	    {
		fieldIndex = 0;

		mass = Double.longBitsToDouble(intBits);
	    }
	    else
	    {
		intensity = Double.longBitsToDouble(intBits);
		pList.addPeak(new Peak(mass, intensity));

		//	System.out.println(mass + " " + intensity);
	    }

	    fieldIndex++;
	}

	return pList;
    }

    public static void  decode64(MzxmlPeakList mpl)
    {
        String input = mpl.getEncodedM2zAndIntensities();
	byte[] byteArr = Base64.decode(input);

	int floatBytes = 64 / 8;
	int fieldIndex = 2;
	int i;
	double mass=-1;
	double intensity=-1;
	PeakList pList = mpl;

	if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	for (i = 0; i < byteArr.length - floatBytes + 8; i += floatBytes)
	{
	    long intBits = 0;
	    intBits |= (((long) byteArr[i]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 1]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 2]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 3]) & 0xff);
	    // Must be in IEEE 754 encoding!
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 4]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 5]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 6]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) byteArr[i + 7]) & 0xff);

	    if (fieldIndex == 2)
	    {
		fieldIndex = 0;

		mass = Double.longBitsToDouble(intBits);
	    }
	    else
	    {
		intensity = Double.longBitsToDouble(intBits);
		pList.addPeak(new Peak(mass, intensity));

		//	System.out.println(mass + " " + intensity);
	    }

	    fieldIndex++;
	}

    }
      
    public static void decode32(MzxmlPeakList mpl)
    {
        String input = mpl.getEncodedM2zAndIntensities();
	byte[] byteArr = Base64.decode(input);

	int floatBytes = 32 / 8;
	int fieldIndex = 2;
	int i;
        if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	float mass=-1;
	float intensity=-1;

	PeakList pList = mpl;
	
	for (i = 0; i < byteArr.length - floatBytes + 4; i += floatBytes)
	{
	    int intBits = 0;
	    intBits |= (((int) byteArr[i]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 1]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 2]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 3]) & 0xff);
	    // Must be in IEEE 754 encoding!

	    if (fieldIndex == 2)
	    {
		fieldIndex = 0;

		mass = Float.intBitsToFloat(intBits);
	    }
	    else
	    {
		intensity = Float.intBitsToFloat(intBits);
		pList.addPeak(new Peak(mass, intensity));

	//	System.out.println(mass + " " + intensity);
	    }

	    fieldIndex++;
	}

    }
    public static void  decode(MzxmlPeakList mpl) throws Exception
    {

        int precision = mpl.getEncodingPrecision();
	if(precision==32)
	    decode32(mpl);
	else if(precision ==64)
	    decode64(mpl);	    
	else
	    throw new Exception("Invalid precision");
	    
    }

    public static PeakList decode(String input, int precision) throws Exception
    {

	if(precision==32)
	    return decode32(input);
	else if(precision ==64)
	    return decode64(input);	    
	else
	    throw new Exception("Invalid precision");
	    
    }
    
    public static PeakList decode32(String input)
    {
	byte[] byteArr = Base64.decode(input);

	int floatBytes = 32 / 8;
	int fieldIndex = 2;
	int i;
        if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	float mass=-1;
	float intensity=-1;

	PeakList pList = new PeakList();
	
	for (i = 0; i < byteArr.length - floatBytes + 4; i += floatBytes)
	{
	    int intBits = 0;
	    intBits |= (((int) byteArr[i]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 1]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 2]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 3]) & 0xff);
	    // Must be in IEEE 754 encoding!

	    if (fieldIndex == 2)
	    {
		fieldIndex = 0;

		mass = Float.intBitsToFloat(intBits);
	    }
	    else
	    {
		intensity = Float.intBitsToFloat(intBits);
		pList.addPeak(new Peak(mass, intensity));

	//	System.out.println(mass + " " + intensity);
	    }

	    fieldIndex++;
	}

	return pList;
    }

    
    public static void decode32(String input, PeakList pList)
    {
	byte[] byteArr = Base64.decode(input);

	int floatBytes = 32 / 8;
	int fieldIndex = 2;
	int i;
        if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	float mass=-1;
	float intensity=-1;

	
	for (i = 0; i < byteArr.length - floatBytes + 4; i += floatBytes)
	{
	    int intBits = 0;
	    intBits |= (((int) byteArr[i]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 1]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 2]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) byteArr[i + 3]) & 0xff);
	    // Must be in IEEE 754 encoding!

	    if (fieldIndex == 2)
	    {
		fieldIndex = 0;

		mass = Float.intBitsToFloat(intBits);
	    }
	    else
	    {
		intensity = Float.intBitsToFloat(intBits);
		pList.addPeak(new Peak(mass, intensity));

	//	System.out.println(mass + " " + intensity);
	    }

	    fieldIndex++;
	}

    }

   
    public static void main(String[] args) throws Exception
    {
	String str = "QzwvlkAuixJDPTc1QNPQ1ENDNKRBB4q5Q0Q/KkDQkWdDRSSFQrfLokNGJhpBgAxQQ0dFuEF7tRNDSBbEQMF19kNJHcRBa3DRQ0swYEFsYMlDTFZUQAAHXkNNGjhAFWy9Q08l3EEqSjRDUB+iP+cgTkNSJfdAQmzOQ1MVKEGQyglDVF53QRGjfENVDexBEDa+Q1X8gUCMB9xDVxuBQW/iQkNYN1RCH3baQ1kuyEFjX+BDWrAwQXXTXkNbYjBBful1Q1wsykBLgTVDXRN4QbyDIUNeRV1AYdtZQ183FkEV2yxDYHBTQY6pPENhQzNBueZAQ2IpDkGUPRNDYyESQfmMDENkHZ9DC1NTQ2UpKkH7iaBDZj0fQKiyNENnJypBpxWuQ2hcVkEPNPVDaiF+Q6GVRUNrHnNCEaYVQ2wVikEClypDbUzwQhqC8ENuQVtBX0QmQ28XUkD8BKZDcC2xQRNMdENyKnRDyaYgQ3M2D0J+S/FDdD22Qh2ZCkN1L/hBrrSSQ3YppUFYJ8RDdzRkQZOySUN4YXZA7hwLQ3k65kEejctDem1qQIXqNEN7RChBy52/Q3wpt0CwN9pDfUv4QOAmvkN+MxxB7V7QQ38vCkCpOI1DgCcIQWa07kOAm1xBbADhQ4EeokC4U71DgaEwQYVtRkOCOxlBpakhQ4KaZUGnmIdDgyAjQVsoJEODqvRBVWZ8Q4QZSkEL5cFDhJljQZtw1UOFn/JBhFz6Q4YhskD4m15DhqKYQcZhIUOHEj5BnH4XQ4ehOEFTJa9DiC5XQcaMWEOItH5BGm9nQ4ksSUEpuyJDiaUWQYvDNkOKI/NA/lviQ4qb6kF5q/hDi0JMQSF8h0OLpFNBsAFyQ4wb80FRQixDjKZMQU5rA0ONFkNBt/coQ42O0EE02fxDjidKQc7AkUOOnUpA9cUYQ48oE0GQwk5Dj5t3QU48t0OQL/JAq2OHQ5CmdkEDogVDkSTZQJGgYEORnPBBQIHCQ5H7Zj/3U4tDknx9Qa8fYEOTGDpAulAAQ5OrgEHAW+ZDlCRrQnBAfEOUqWJBluPJQ5UzdEHdeXlDlZWfQhN3QEOWGDhBwdKkQ5aQIkEs7MtDlx/eQkuVd0OXi3dBM/NrQ5gg8EFV9LxDmKNEQXAPQkOZITVA9iTQQ5mrO0AvOmVDmjukQYPu5kOauYhBMnnUQ5tCSkEHGUpDm5yOQUNT50OcmQBEFHJ2Q50agkKjIudDnawTQbYpokOeHexBGf+dQ56yZEKROY5Dnx9CQW++dEOfqCJBg+u2Q6AzxEGjs3xDoKt0QZPgjkOhKk5BTM22Q6GnhEF24mVDoiiJQMeQREOioSRCnBpMQ6MZ7EFDQHxDo6YKQZbj6kOkKtRB/n7uQ6ScFkG2BDlDpSqKQcWjAEOlovhB04jcQ6YafEEKSplDpqkuQbcPdEOnK5JA+E9YQ6ekZEGJLapDp/7TQNmRW0Oomf9CDxn1Q6kmQkKU7qhDqa5/QoHI70Oql89EBAlkQ6sgnkK6CbJDq6pYQYAX1kOsRfZBkXW0Q6zAPkGODatDrSbaQYfvC0Ottj9CvpfcQ648dkCKAgVDrqHuQXz5wUOvCFY//3hoQ6+0gEIFvzJDsD64QklsZEOwlXBDX1ZRQ7EhukJ84ZRDsb6/QdWLlUOyKadCUOXxQ7KZQUHjK75DsxeRQi/UX0Ozp85CkYkXQ7Qb9kEGZMZDtKwsQdeqCkO1GrJAP8qHQ7WhwkG2uv5DthtOQQ0zt0O2pKJCAKWKQ7crSkDukwRDt6EWQcLR3kO4OoZBN3ewQ7ixxEITLx9DuS6eQoEJ90O5ou5B1aV8Q7pNxkHxWvFDuqRMQkqgrkO7HmBBr0hmQ7umkEIc8LhDvFGKQSILeEO8xEJBR4ENQ70bjEGqedtDvaDQQg6BdkO+IXBBPhz+Q76gNkIWoSdDvwnzQeWpvEO/m3BAudtnQ8AieUFIYPhDwJ4OQkisKkPBRIZBSfnkQ8GqHEHppm9Dwi3mQZkQ5UPCqehBnZWkQ8M0/EEBqbBDw7SsQQc47kPEaWRBSTRCQ8T5J0F+mq1DxVAjQScjyEPFn/5CBOc3Q8Y10kGXIUxDxq27QjPxHkPHJ9FBeJmWQ8ejvEIjHBJDyC9nQWwlykPIu4hBrDQKQ8k1d0GaVTRDyagGQjIUKEPKJ0NBD2N3Q8rS1kEAER5DyyruQRhtuEPLnwhCVtCnQ8wlQkIQdWxDzK7xQjG6bUPNOgtCFlAHQ82v2EGBPldDzkRGQk7OHkPOn/pCKaAJQ88kakHuhK9Dz6HIQjymo0PQObhCHx8QQ9CzgkIMrCND0T08QbT2WEPRr1BBYqqSQ9IYvkCHQBVD0stmQZ/D0kPTSCZBL+1AQ9O620IN4LBD1CsUQyf28EPUoapCgrt0Q9VE3EKj06JD1ad/QrC86kPWJ/VCMwKfQ9aevkJESwFD1xojQldow0PXsdxB6hvTQ9hNg0FaeLtD2MaVQiKQj0PZPiJDBq/0Q9m/yEH85ktD2i5RQSvdyUPasd1CGl/5Q9s0FEESGTND26xMQYJihkPcoItC5sv4Q90la0IflAtD3ar8QTZ9PUPeC9FCifhYQ96zjUJRLoND3znDQk4iPkPfq81B2AYPQ+AzHUGgxMBD4IdPQUnF40Pg1kBBuDuRQ+FBOEGHSx5D4bnTQv0h1kPiKU1DUZoXQ+KuB0LQFs1D4y2aQkbWPUPjo+JCGr2IQ+QeokJPnk5D5LBNQbX870PlHyRBmA70Q+Weq0IZ+cdD5jSPQdDqrEPmr7lBUrswQ+ceiEGGo/pD56MOQZ6ZnkPoNApB30d/Q+i1wUIzcqJD6TdMQd1UOUPpoZNBFIgDQ+paIEJggFhD6ql8RAPnRkPrOKxCqB6ZQ+vOaUHxz29D7CaXQxQEp0Pssy1CQDJoQ+035kGczalD7afoQWx8NEPuMgpA7xvdQ+7FW0E0hypD731zQ3SC/EPwDKRCAUo7Q/CnKUI1CExD8Tc7QgAhl0Pxw/9CZvyQQ/IgIEOwPTBD8qnQQuRKFEPzJdhCtyyQQ/PC10MQqaRD9DJUQnmWjUP0qeVCWgWWQ/VJ8UJpPuRD9bXJQy4T9UP2L8NCNgFMQ/bYMEIk9fdD91hoQYzOGEP3vdhBYV1oQ/hA9kGG0FRD+LeXQdhqcUP5LxBCUvvAQ/mnO0JkEzVD+jTeQlUZw0P6sQBCb4IMQ/sTDUHvIJZD+8LDQYSrj0P8LY5B/MqEQ/yfFkHIeSJD/STTQdXzukP9we5DjR+KQ/5Ae0IMcQJD/rT7QjaN6UP/TnZBvBBsQ/+qf0H4w21EAAj+QgOwAEQAVsRDSri9RACYI0KCMHVEAM9aQkAAVkQBCvBCKxXrRAFPpEJFB6REAYktQdksE0QB5wpBzrsZRAIqKEHEizhEAlp6QdonEEQCjyZCNZQBRALXUUI7AAhEAx4nQqgrFEQDVzpC1QPcRAORkkMcRNJEA9InQowkW0QEHExCVJY2RARW00I6lyREBNg3RZ9ws0QFGXdEleh5RAVankMIEZJEBZmYQhuRaUQF0zZCRB9aRAYJuEI1v3pEBmMmQ4ef2EQGny5CNpsIRAb0YkOKUHpEByz6Q3OuIUQHYyJCTj0ZRAepL0LWUzREB9rsQt17vkQIGnpCjn+YRAhatkJd5CRECJMUQb1ALkQJEOREEYVVRAlOWESMlM9ECYWkQ5fH6EQJy3FCiU69RAoVSEI0dblEClJYQiMfR0QKq69Cmso+RArZJ0NMk/tECzPgQt5IckQLaTRDc5uYRAubx0Kd7sxEC8bgQi6SX0QL8V5BTPS2RAwd0EHRXNNEDF2AQPAqUkQMpL5CBU3NRAzkHkJ7XXJEDSSWQslhiUQNc5tB/7rKRA2cUUI9YOFEDd8uQoeiQkQOKvZBoYHzRA5g5EISL69EDpbeQh9V/0QO2uRCNfDORA8SiEME3blED1RCQuI4OkQPlN5BNX9mRA/W40E32GpEEBB8Qez/EEQQXv9B5OeuRBCOwEHBvKdEENprQlKFtkQRFE9BaBQ2RBFgekJ+aplEEaDoQzQXJkQR1fFDjza/RBIZM0N4J9tEElsqQ5dm90QSovRCbZ9YRBLVhkIBSvFEExQOQkSSokQTWjRCJERkRBOiOkFg1o1EE9cCQiarlEQUK4xC1AwsRBRdWUIoxX9EFKqZQm3pX0QU235CVRlLRBUhekLCEMlEFVNcQ9+9LkQVkp9DINZMRBXmJEL5YWpEFhXwRFQz10QWVj5Dl7NzRBadGEWFap5EFtvQRJWqdUQXHcRCyCnqRBdtikJuq3hEF58cQmsCe0QX01hBv/LyRBga6kJmahZEGF3HQmqFDkQYnbJCB/Z4RBjSc0Jau1NEGRuAQs8K90QZXuZDXyQhRBmXIUOAnHVEGdeOQzDmA0QaFQJB5dBARBpKwkI6xd1EGn7xQZ51SEQaq7BB9y/eRBrcsUG+BFREGxuMQjI/IkQbWeRCcPYURBumPkJoJA9EG95oQuMRPUQcHgdCRW2xRBxg4kL+UD9EHKUGQq9V20Qc2ZNCdZzgRB0NbkMJuLVEHV5ZQmtdh0QdpbFCn4wgRB3VkEMTGetEHhUYQmZtS0QeUyBCCaTpRB56gEG331REHqmGQqMDVkQe6lhC6UTbRB8Y9kHSP6lEH1CIQo+wvkQfnzZDmQpQRB/ZSEO/6xxEIBuCRDtZaUQgXyJDqCpJRCCctkMf4+ZEIOb8QkkzWEQhK3ZCjyguRCFjwkINUJNEIZWmQbFfHkQh59JDkfJCRCIvVEROlZ1EImEOQysnhUQio6RBAyVORCNE5kCrShREJT2jQFOOG0QmX+ZBjbR8RCaVpEDf38xEJxRfQ6QvlkQnVItCu0gvRCfUvkNbJZpEKBdOQrUcPEQoX7RCZqMgRCimgEGQuXlEKOHyQa/pBEQpF5ZAlGaQRClnxEECTpxEKaIQQYOQuEQp60pBJ66sRCobS0LFS1lEKlQkQk8eEEQqnSpCe262RCraWEPbyvZEKyN8Qxhp4kQrXOZDd+kLRCuYKEMBlTBEK9uWQZTzCkQsE8hAATgVRCxYbEKFwvlELJpEQdaU7kQs20JBhY85RC0SCkHMRIVELV5WQZYbYkQtkDRBmMEARC3fgkNmaNpELht8Q6NZe0QuVa5Cum7KRC6VjkIptZZELti2Q0gwgEQvGVBCvPcFRC9d6EIosPlEL6EmQeplMUQv5EpB4FFjRDAf7EIC60VEMGq+QYK7JkQwqX5CL1o+RDDjLEGkHn9EMRNgQgdKmkQxYvBBieTARDGd5EDqk61EMeCOQbFVNUQyWLZE3MqfRDKbBkPyuZBEMtzsRRk+p0QzHoRES3CSRDNkfELTQHBEM5geQMAEg0Qz35RBoZJDRDQqYEC8ZipENFe+QPE5skQ0nQxBp7VhRDTaDEGpJ5pENSd6QQv7HkQ1T5ZBOKCrRDWmDkFmuEBENdzuQVPdyEQ2E7JBkZ/GRDZ1tEEQ7AZENs0cQX7BGEQ3ICJBBByARDdYjkDeplNEN6LEQNuHREQ33YJBNn6FRDgYWEEZgRhEOFwGQXa8nkQ40VxDGvy+RDkTqkKcM5pEOVkmQlZgEEQ5oSBBgfy+RDnjcEIArBJEOiDMQUcsZEQ6W6pBTMgQRDqiwkFo93JEOtvuQZGgPkQ7Ik5CDIykRDtYrkJRkaFEO3+0QGnqykQ7ptBBZLurRDvsbkCaKOlEPD9EQZ7fGkQ8eHpBZdJoRDyzbkES+PpEPQ6yQzpD80Q9WHRCq98yRD2iUkJU23ZEPec8QgvO7EQ+G+xCMRtbRD5cykKLZEpEPqB6QZTJnEQ+3BRBjHdKRD+cMkOj8thEP9fcQ6vkSURAGyRC10EdREBeGEIAF1NEQJs6QUhp50RAzEZBsdSvREErqkGBgMREQVuaQWVWxkRBmVpBhTY7REHcXkD5+UNEQij8QXnw+ERCXRxCaW/sREKnakHPj81EQtjkQgZI+ERDEJxCCyuYRENdkEExfbBEQ5HCQYecrEREFWJFHkAfRERXNkRU2lRERJnyQy/+cURE5nJB5O0HREUxHkGC7qdERWKOQaj6c0RFlBBApMXJREXk+kD7WchERh6kQQMqNERGW7ZCsRJXREagSkIsWntERvEMQb0XakRHLTJBr3e3REdkVEGgbHpER6V+QgFx0kRH6JpB0I0EREgmEkHI4RlESGFQQedOhkRInw5BeDXIREjWqkFezelESR5+QUCKGkRJYWpAT4IuREmfpkITa0NESdz6QowaGkRKNlRAkNK0REpdhEHJ4pdESp6EQcxITkRK4BpBVxsWREsisEHvRBBES1/UQYlbN0RLmCZA4MKrREvYnkLJuuxETBaQQhK0D0RMT+hBbHVYREyC/EGWxIpETNkWQQXE4ERNEMZBBujrRE1PREGH1eFETZi4QdG9gkRN41RBj+EPRE4olkGrSX9ETlaIQr8Um0ROnmBCKa4+RE7vCkBN3hNETy9AQooBiERPYN5CAT+eRE+dTkGVxklET8tSQMEUYkRQGCRBOqQ/RFBHwkCLgtBEUKeaQWtvgkRQ7kpBzWfWRFEgrkDNTFFEUXAwQVlp/URRn6pAighRRFHefEFtiOxEUiJWQQ5Lk0RSXd5AwfUdRFKhHEGPxdNEUtu8QXJvuERTGURBO2egRFOdVkUFd/VEU98WRIBXQERUIJZDVVvSRFRkdkGSHUFEVKfAQZrqKURU285BwBpURFUNhkHo/bZEVVncQkO2A0RVmkhB498cRFXbHEGOC8pEViVQQJRqrkRWWExBFA3LRFaWZEEZpvZEVtIIQScGKERXH9RBNLVWRFddiELERkhEV58QQoX0mURX2UBCioFZRFgbKkKH+hBEWFzEQhck1ERYngRA9v0kRFjOYkDFwzxEWR7AQqqUwERZX3JCR1koRFmPmkHC+dhEWdB6QdX9oERaI0RBnFo0RFp+EEE++L5EWsZEQY7EfkRbFcxAupVeRFtOpEKVIp1EW3qCQlwDokRbp5pCUnb5RFvn+EJBLK5EXFX8QqrWUkRclkJBOFzNRFzSyEEvRbREXRYqQaYC60RdXpxCY2RiRF2bukIcosZEXd6AQa0ke0ReFkBBQwx0RF5VOkC4oSJEXo5GQRriwkRe7f5Bs250RF80skEn/9BEX3B4QWrSLkRfpIBBKr5URF/aCkDJtH1EYDAmQTDXhERgZrhBUmFKRGChKEGBOPlEYMvEQTvbHERhF2RCTzjSRGFGbkIs0B9EYa4uQKuSxkRh125BRsEMRGIhXECa1W9EYmE8P/RlXERim8pBHWeuRGLfYkIW9UFEYyEGQmqHx0RjZPpCZgiCRGOnHEFKxxBEY+b8QbsDgURkIKRBbsWRRGR48kIWIrZEZLL0Qd1iiERk7spBqmKCRGUb+kHnBUFEZWPqQfAtr0Rlp/RBjbm5RGXT5ECHVSpEZhwUQBhz4kRmeeBBjIU6RGaxZkEM1XxEZuAwQUtsTkRnFl5BGJywRGdgzEIGkVVEZ5s6QcdnjERn3kZBlBqZRGgzlkFbNNBEaJRCQcWZSkRo2f5AFqOyRGkUQEE7KbxEaXg0QN+xEURptpZBCbGzRGnjOEDObtFEaiiUQZxClkRqZx5BYp0YRGqf6kDyPBpEatOEQTfd6URrGZ5BgbmORGtV1EFP8kVEa5nkQURcr0Rr2wxBE1rMRGwxAEEPsHZEbGWiQRVOqERs2opBuasSRG0jaEGjMI1EbWCkQkt1a0RtmLxCKRpzRG3U9kIQRC5EbfvCP/Eq2ERuLnxBI4CORG5opkFVmapEbpZkQPM96kRu2a5Bg0WfRG8q1EMjK/1Eb2A0Q8NP8ERvo4pDKIMGRG/goEJJHCREcBekQbe+EURwaBhBA2WJRHCspEC/RChEcNvaQN8GIURxN2ZBjPZIRHFkREHTc8NEcaKOQVDgKERx25JCFZzHRHIgzkLWDVhEcmBwQgmpjkRypLRA2yNYRHLhjEHiwSdEc58ORQ7Ug0Rz3vRElBayRHQhLkN6eYlEdFd+QYo2NER0lsJBxfUURHTb3kGi4BVEdR+8QRadCER1cYBAzVVMRHXlakH/1jJEdiSsQoQ++0R2YWRBlEHwRHaTkkGY3vBEduHaQQddZ0R3JFRBiO1/RHdVskFc+IREd5SwQSDV+0R4Fd5BUTC8RHhdAkFdnBNEeKHEQaWoBUR44WxCtRtXRHkZhEK9YnFEeWLMQpw+xER5n7RBnAJbRHneWEAn1/ZEeiiqQaTxWER6a/xAFq9hRHr++EAzyohEe6YaQOanGkR71ZZAKahhRHxM7EBsalVEfMdwQKCjhER9HeBBtNUSRH1fjEMYj8hEfabEQpZoL0R95sBAifW0RH5nkEGaepVEfqAaQDbtRER+33JBLriMRH9cEEGE9FBEf60+QPWqN0R/3nBAvIbpRIAS6kB1P9tEgDJLQEkfxkSAWMJA+vO1RIBvx0CjCk5EgJWQQUc2EESAsC1B6+uJRIDQyUGEEsBEgPF0QBmqu0SBMD1BMe2+RIFWjkFQT0dEgW6iQXntgESBkUBA5r/9RIGxvkDgKCREgdC6QYA4nkSB6aJBUaOHRIH82UELQMNEghGtQR8B90SCOyRBeLMPRIJca0GDy6xEgnm3QKYthUSCkJZATRcuRIK08EDU37ZEgtWeQakzJkSC8gpBgMz4RIMNIUFf6UBEgy+tQRg8ckSDVdpAu5TyRIN0ykFuKEpEg49aQLhg8USDsu9CD6xaRIPOzEDIpixEg+e9QHpSikSEFeZBJsonRIQrX0D6DUVEhFULQBehuUSEaPpBmJK/RISMSkFdp6lEhK3vQSKh80SE7sVDmPDWRIURwkMV9lZEhSx8QqKgmUSFTlVCBwQmRIVzd0HaB/ZEhY6hQF4UO0SFrlNBQZhoRIXNOkCXi31Ehe5eQMvAkkSGMCVDVwekRIZSuUL5cxpEhnbVQh2nTkSGlehAorzhRIat5kFfw9JEhs0KQYyyFUSG63lBxAmZRIcuUEQN1QpEh07BQ47+zkSHbytCx4ijRIeb+EFdZBREh7HOQUM0bUSH1cVBNdc7RIf1RkID0fJEiBVtQVJ0dESINCtADvvURIhu/0BxcTxEiI5JP+bcK0SIrwhAVFCVRIjMGUCNgPBEiPuIQc7yfESJD9ZCmS9oRIktLUKRiiVEiU7BQZdup0SJafVA8NQ2RImK4kDJMIJEia4IQKKQMkSJ8j5BFQSqRIoTr0CQ6T5EijBbQJYoxUSKbfRAlN7KRIrhE0EvKV5EivfiQRUQC0SLGTJBCJUcRIs67UGJeXNEi1fGQZUJxESLdehB+ca3RIudEEEoaahEi7U9QJ2gXESL041AxqJVRIvw0EDwVHlEjCnWQRjy/USMTERBrjtERIxkj0FmfdpEjHlnQCI7WESMmk1BGJ5MRIyz4EEhpDhEjNftQKcDokSNEXNBGtwwRI0sHT/o6yREjZCjQOgrcUSNt2xBH5ayRI3PVkED6k9EjhvYQHxM/ESOMDBAvnUgRI5MS0GDPqREjnR4QacuTkSOkWhBzDvRRI6wBUHnae9EjtCwQacszESPEoJCI63VRI8vg0CH6k5Ej3QOQI7uFkSPj4VAgMifRI/KiUJsux9Ej/OtQerjOkSQEcFA91V8RJBEPEBdDFdEkKzFQFbA/ESQz3hAjEu4RJD0rUAwEr5EkUdzQChm1kSRWy1AYSolRJFwK0C8d/JEkZs/QTMeDESRuS1BjmUDRJHXVUDdY15Ekfw4QgVQJ0SSEP1DKFWkRJIxykKT58NEklDLQVnmP0SSjjJAVRahRJNQ2kCIV4tEk+GwQJGbNkSUEdhBNmyaRJQ3VEC30zZElFE4QFwJoUSUdB9ARO7dRJS5UECKeu9ElO/IQIWA+USVNppAnD7bRJVKOkEPpW9ElXVOQCjqAkSVvP5Ac+pMRJXQ70B5sl5Elhs8QJRoA0SWNiNCGkTrRJZXvUHV4sZElnehQYrbpkSXUWtABT7zRJdwLD/bU0ZEmYUJQKBLs0ScuwNAol8URJ0SyUAZWIlEnZJPQA2nEkSeNZZAhfhaRJ7QMUCuZ+tEn4ujQJHlyUSfzi9AX5ZERKJ/YkBo5B9Eo14jQGm3mESjdmhANGpdRKPWi0ByaUJEpBRDP/Ma5kSouJFAau/sRKkPMUAuDlFEwDB5QBz8GUTI1aBACewbAA==";
	PeakList pList = MZXmlHandler.decode(str, 32);


    }
    
//    public static void decode(byte[] input, int precision)

}
