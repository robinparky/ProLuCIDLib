

/**
 * @file ScoreCalculator.java
 * This is the source file for edu.scripps.pms.mspid.ScoreCalculator
 * @author Tao Xu
 * @date $Date: 2005/10/11 20:51:50 $
 */
package edu.scripps.pms.mspid;

public class ScoreCalculator {


    /*
     * hypergeometry calculates the hypergeometrical
     * probability: numTotal elements, numMarked of them marked.
     * numChosen elements are drawn from the pool.
     * what it is the probability that numMarkedChosen elements
     * are marked.
     * numTotal >= numChosen >= numChosen-numMarkedChosen; numTotal >= numMarked >= numMarkedChosen;
     */
/*

    public static float hypergeometry(long numTotal, int numMarked, int numChosen, int numMarkedChosen) {
System.out.println("in hypergeomtry, numTotal: " + numTotal + "\tnumMarked: " + numMarked + "\tnumChosen: " + numChosen + "\tnumMarkedChosen: " + numMarkedChosen);
        double ftemp = 1.0;
        int numUnmarkedChosen = numChosen - numMarkedChosen; // num not marked in the draw
        long numUnmarked = numTotal - numMarked; // num not marked in the pool
        int numMarkedNotChosen = numMarked - numMarkedChosen; // num marked not chosen
        long numNotChosen = numTotal - numChosen; // num in the pool not chosen
        if(numNotChosen < 0 || numMarkedNotChosen < 0) {
            System.out.println("WRONG input to HYPERGEOMETRY");
            return 0.f;
        }
        for(int i = 1; i <= numChosen; i++) {
            if(i <= numMarkedChosen && i <=numUnmarkedChosen) {
                ftemp *= (numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen+i);
                ftemp *= (numMarkedNotChosen+i)/(double)i;
            } else if (i <= numMarkedChosen && i > numUnmarkedChosen) {
                ftemp *= (numMarkedNotChosen + i)/(double)(numNotChosen + i);
            } else if (i > numMarkedChosen && i <= numUnmarkedChosen) {
                ftemp *= (numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen + i);
            } else if (i > numMarkedChosen && i > numUnmarkedChosen) {
                ftemp *= (double)(numNotChosen + i);
            } else {
                System.out.println("WRONG Configuration");
                System.exit (1);
            }
        }
        if(ftemp < 0.0 || ftemp >= Float.MAX_VALUE) {
            double dtemp2 = 1.0;
            ftemp = 0.0;
            for(int i = 1; i <= numChosen; i++) {
                if(i <= numMarkedChosen && i <=numUnmarkedChosen) {
                    dtemp2 = 1.0;
                    dtemp2 *= (numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen+i);
                    dtemp2 *= (numMarkedNotChosen+i)/(double)i;
                    ftemp += Math.log(dtemp2);
                } else if (i <= numMarkedChosen && i > numUnmarkedChosen) {
                    ftemp += Math.log((numMarkedNotChosen + i)/(double)(numNotChosen + i));
                } else if (i > numMarkedChosen && i <= numUnmarkedChosen) {
                    ftemp += Math.log((numUnmarked-numUnmarkedChosen+i)/(double)(numNotChosen + i));
                } else if (i > numMarkedChosen && i > numUnmarkedChosen) {
                    ftemp += Math.log(i/(double)(numNotChosen + i));
                }
            }
            ftemp = Math.exp(ftemp);
        } else {
            return (float)ftemp; 
        }
        if(ftemp < 0.0 || ftemp >= Float.MAX_VALUE) {
            System.out.println("Failed Accuracy");
        }
        return (float)ftemp;
    }
*/
    /*
     * HYPERGEOMETRY calculates the hypergeometrical
     * probability: numTotal  elements, numMarked of them marked.
     * numChosen elements are drawn from the pool.
     * what it is the probability that numMarkedChosen elements
     * are marked.
     * numTotal  >= numChosen >= numChosen-numMarkedChosen; numTotal  >= numMarked >= numMarkedChosen;
     */
    public static float hypergeometry(long numTotal, int numMarked, int numChosen, int numMarkedChosen) {
//System.out.println("in hypergeometry, numTotal : " + numTotal  + "\tK: " + numMarked + "\tnumChosen: " + numChosen + "\tnumMarkedChosen: " + numMarkedChosen);
	double ftemp = 1.0;
	int J = numChosen - numMarkedChosen; // num not marked in the draw
        long L = numTotal - numMarked; // num not marked in the pool
        int dK = numMarked - numMarkedChosen; // num marked not chosen
        long dN = numTotal - numChosen; // num in the pool not chosen
	if(dN < 0 || dK < 0) {
	    System.out.println("WRONG input to HYPERGEOMETRY");
	    System.out.println("N = " + numTotal +"\t"+"K= " + numMarked+
			       "\t"+"numChosen= " + numChosen + "\t"+"numMarkedChosen = " +numMarkedChosen);
	    return 0.f;
	    //System.exit (1);
	}

	int i; // counter
	for(i=1; i <= numChosen; i++) {
	    if(i <= numMarkedChosen && i <=J) {
		//ftemp *= (double)((dK+i)*(L-J+i))/(double)
		//  ((dN+i)*i);
		ftemp *= (double)(L-J+i)/(double)(dN+i);
		ftemp *= (double)(dK+i)/(double)i;
	    }
	    else if (i <= numMarkedChosen && i > J) {
		ftemp *= (double)(dK + i)/(double)(dN + i);
	    }
	    else if (i > numMarkedChosen && i <= J) {
		ftemp *= (double)(L-J+i)/(double)(dN + i);
	    }
	    else if (i > numMarkedChosen && i > J) {
		ftemp *= (double)i/(double)(dN + i);
	    }
	    else {
		System.out.println("WRONG Configuration");
		System.exit (1);
	    }
	}
	//System.out.println("ftemp = " + ftemp);
	/*
	 * If the value is outside of accuracy try to recalculate it
	 *
	 */
	if(ftemp < 0.0 || 
	   ftemp >= Float.MAX_VALUE) {
	    double dtemp2 = 1.0;
	    ftemp = 0.0;
	    for(i=1; i <= numChosen; i++) {
		if(i <= numMarkedChosen && i <=J) {
		    dtemp2 = 1.0;
		    dtemp2 *= (double)(L-J+i)/(double)(dN+i);
		    dtemp2 *= (double)(dK+i)/(double)i;
		    ftemp += Math.log(dtemp2);
		}
		else if (i <= numMarkedChosen && i > J) {
		    ftemp += Math.log((double)(dK + i)/(double)(dN + i));
		}
		else if (i > numMarkedChosen && i <= J) {
		    ftemp += Math.log((double)(L-J+i)/(double)(dN + i));
		}
		else if (i > numMarkedChosen && i > J) {
		    ftemp += Math.log((double)i/(double)(dN + i));
		}
	    }
	    ftemp = Math.exp(ftemp);
	}
	else
	    return (float)ftemp;
	
	if(ftemp < 0.0 || ftemp >= Float.MAX_VALUE) 
	    System.out.println("Failed Accuracy");

	return (float)ftemp;
    }
    public static double nChooseK(int n, int k) {
        double result = 1;
        int nMinusK = n - k;
        int smaller = k;
        int greater = nMinusK;
        if(k > nMinusK) {
            greater = k;
            smaller = nMinusK;
        } 
//System.out.print("n: " + n + "\tk: " + k);
        while(smaller > 0 ) {
            result /= smaller--;
            if(n > greater) {
                result *= n--;
            }
        }
//System.out.print("\t" + result);
        while (n > greater) { 
            result *= n--;
        }
//System.out.println("\t" + result);
        return result;
    }
    public static void main(String args[]) {
        for(int i= 1; i < 1000; i++) {
            for(int j = 0; j <= i; j++) {
                nChooseK(i, j);
            }
        } 
    }
}



