/**
 * @file SpectrumReader.java
 * This is the source file for edu.scripps.pms.util.stats.StatCalc class
 * @author Tao Xu
 * @date $Date
 */



package edu.scripps.pms.util.stat;

     /* 
        An object of class StatCalc can be used to compute several simple statistics
        for a set of numbers.  Numbers are entered into the dataset using
        the enter(double) method.  Methods are provided to return the following
        statistics for the set of numbers that have been entered: The number
        of items, the sum of the items, the average, the standard deviation,
        the maximum, and the minimum.
     */
     
     public class StatCalc {
     
        private int count;   // Number of numbers that have been entered.
        private double sum;  // The sum of all the items that have been entered.
        private double squareSum;  // The sum of the squares of all the items.
        private double max = Double.NEGATIVE_INFINITY;  // Largest item seen.
        private double min = Double.POSITIVE_INFINITY;  // Smallest item seen.

        public static void main(String[] args) {
            double z = 1;
            z = 1;
            System.out.println(z + "\t" + zScore2PValue(z));
            z = 1.65;
            System.out.println(z + "\t" + zScore2PValue(z));
            z = 1.96;
            System.out.println(z + "\t" + zScore2PValue(z));
            z = 2.58;
            System.out.println(z + "\t" + zScore2PValue(z));
            z = 3;
            System.out.println(z + "\t" + zScore2PValue(z));
            
            double t = 0; 
            double df = 0;
            t = 5;
            df = 1;
            System.out.println(t + "\t" + df + "\t" + getTTestPValue(t, df));
            t = 5;
            df = 2;
            System.out.println(t + "\t" + df + "\t" + getTTestPValue(t, df));
            t = 2;
            df = 3;
            System.out.println(t + "\t" + df + "\t" + getTTestPValue(t, df));
            t = 2;
            df = 20;
            System.out.println(t + "\t" + df + "\t" + getTTestPValue(t, df));
            t = 2;
            df = 1000;
            System.out.println(t + "\t" + df + "\t" + getTTestPValue(t, df));

	    // assume all the values are greater than or equals to 0
	    double[] darr = new double[5];
	    darr[0] = 1.0;
	    darr[1] = 2.0;
	    darr[2] = 3.0;
	    darr[3] = 1.5;
	    darr[4] = 0.4;

//	    int[] arr = getHistogram(darr, 1, 10);
//	    for(int i=0;i<arr.length;i++)
//		System.out.println("-->>" + arr[i]);
	   
	    StatCalc cal = new StatCalc(darr);
		System.out.println("-->>" + cal.getStandardDeviation());

            for(int i = 1000; i < 100000; i++) {
                double z1 = i/1000.0;
                double p = zScore2PValue(z1);
                if(p > 0.0499 && p < 0.0501) {
                    System.out.println("Z: " + z1 + "\tp: " + p);
                }
                if(p > 0.000099 && p < 0.000101) {
                    System.out.println("Z: " + z1 + "\tp: " + p);
                }
            }
	    
        } 
        public StatCalc(double [] values) {
            for(double d : values) {
                enter(d);
            }
        }     
      
        public StatCalc() {}
        public void enter(double [] values) {
            for(double d : values) {
                enter(d);
            }
        }     
        // covert z score of a normal distribution to tailed p-value
        public static double zScore2PValue(double z) {
            // Returns the two-tailed standard normal prob. level given a z.
            double absz = Math.abs(z);
            double a1 = 0.0000053830, a2 = 0.0000488906, a3 = 0.0000380036,
	           a4 = 0.0032776263, a5 = 0.0211410061, a6 = 0.0498673470;
            double p = (((((a1*absz+a2)*absz+a3)*absz+a4)*absz+a5)*absz+a6)*absz+1;
            p = Math.pow(p, -16);
            return p;
        }
        
        public static double getTTestPValue(double t, double df) {
            double abst = Math.abs(t);
            double tsq = t*t;
            double p;
	    //var abst = abs(t), tsq = t*t, p;
            if(df == 1)
                p = 1 - 2* Math.atan(abst)/ Math.PI;
            else if(df == 2)
                p = 1 - abst/ Math.sqrt(tsq + 2);
            else if(df == 3) 
                p = 1 - 2*(Math.atan(abst/ Math.sqrt(3)) + abst* Math.sqrt(3)/(tsq + 3))/ Math.PI;
            else if(df == 4) 
                p = 1 - abst*(1 + 2/(tsq + 4))/ Math.sqrt(tsq + 4);
            else {
                // finds the z equivalent of t and df st they yield same probs.
                double z = t2z(abst, df);
                if (df>4) 
                    p = zScore2PValue(z);
                else 
                    p = zScore2PValue(z); // small non-integer df
            }

            return p;
        }
        public static double t2z(double t, double df) {
            // Converts a t value to an approximate z value w.r.t the given df
            // s.t. std.norm.(z) = t(z, df) at the two-tail probability level.

            double a9 = df - 0.5;
            double b9 = 48*a9*a9;
            double t9 = t*t/df;
            double z8;

            if (t9 >= 0.04) 
                z8 = a9* Math.log(1+t9);
            else  
                z8 = a9*(((1 - t9*0.75)*t9/3 - 0.5)*t9 + 1)*t9;
            double p7 = ((0.4*z8 + 3.3)*z8 + 24)*z8 + 85.5;
            double b7 = 0.8* Math.pow(z8, 2) + 100 + b9;
            double z = (1 + (-p7/b7 + z8 + 3)/b9)* Math.sqrt(z8);

            return z;
        }

        // assume all the values are greater than or equals to 0
        public static int [] getHistogram(double [] values, double intval, int max) {
            int num = (int)(max/intval);
            int [] counts = new int[num+1]; 
            double adj = 1/intval;
            max = (int)(max*adj);
//System.out.println(intval + "\tmax: " + max + "\tnum: " + num);
            for(double value : values) {
                int index = (int)(value*adj+0.5);
                index = index < max? index : max;
                index = index < 0? 0 : index;
                counts[index]++;
            } 
            return counts;
        }
        public void enter(double num) {
              // Add the number to the dataset.
           count++;
           sum += num;
           squareSum += num*num;
           if (num > max)
              max = num;
           if (num < min)
              min = num;
        }
     
        public int getCount() {   
              // Return number of items that have been entered.
           return count;
        }
     
        public double getSum() {
              // Return the sum of all the items that have been entered.
           return sum;
        }
     
        public double getMean() {
              // Return average of all the items that have been entered.
              // Value is Double.NaN if count == 0.
           return sum / count;  
        }
        public static double getStandardDeviation(double [] values) {
            int num = values.length;
            double sum = 0;
            for(double d : values) {
                sum += d;
            }
            double mean = sum/num;
            double sumDiffSqr = 0;
            for(double d : values) {
                double diff = d - mean;
                sumDiffSqr += diff*diff;
            }
            return Math.sqrt(sumDiffSqr/(num-1));
        } 
       
	//get standard deviation from frequency distribution
	public static double getStandardDeviationFromDist(double[] values, double mean) {

	    
	    return -1;
	}
	
        public double getStandardDeviation() {  
             // Return standard deviation of all the items that have been entered.
             // Value will be Double.NaN if count == 0.
           double mean = getMean();
           return Math.sqrt( squareSum/count - mean*mean );
        }
        
        public double getMin() {
             // Return the smallest item that has been entered.
             // Value will be infinity if no items have been entered.
           return min;
        }
        
        public double getMax() {
             // Return the largest item that has been entered.
             // Value will be -infinity if no items have been entered.
           return max;
        }


     
     }  // end class StatCalc



