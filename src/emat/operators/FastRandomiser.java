package emat.operators;

import beast.base.core.Log;
import beast.base.util.Randomizer;

public class FastRandomiser {

    
    public static int randomChoicePDF(double[] pdf) {
        double U = nextDouble() * getTotal(pdf);
        for (int i = 0; i < pdf.length; i++) {

            U -= pdf[i];
            if (U < 0.0) {
                return i;
            }

        }
        for (int i = 0; i < pdf.length; i++) {
            System.err.println(i + "\t" + pdf[i]);
        }
        throw new Error("randomChoiceUnnormalized falls through -- negative components in input distribution?");
    }

    
    private static double getTotal(double[] pdf) {
    	double sum = 0;
    	for (double d : pdf) {
    		sum += d;
    	}
    	return sum;
	}


	// essentials from java.util.Random without synchronisation
    private static final double DOUBLE_UNIT = 0x1.0p-53; // 1.0 / (1L << 53)
    private static final long multiplier = 0x5DEECE66DL;
    private static final long addend = 0xBL;
    private static final long mask = (1L << 48) - 1;
    private static long seed = Randomizer.getSeed();
    
    public static double nextDouble() {
        return (((long)(next(26)) << 27) + next(27)) * DOUBLE_UNIT;
    }

    public static int nextInt(int bound) {
        int r = next(31);
        int m = bound - 1;
        if ((bound & m) == 0)  // i.e., bound is a power of 2
            r = (int)((bound * (long)r) >> 31);
        else {
        	// compensate for somewhat non-randomness of r
        	// see @java.util.Random.nextInt for details
            for (int u = r;
                 u - (r = u % bound) + m < 0;
                 u = next(31))
                ;
        }
        return r;
    }

    private static int next(int bits) {
        seed = (seed * multiplier + addend) & mask;
        return (int)(seed >>> (48 - bits));
    }

    
//	static double nextDouble() {
//    	// return FastRandomiser.nextDouble();
//        return wyhash64() * 0x1.0p-53;
//    }
//    
//    static long wyhash64_x = Randomizer.getSeed();//System.currentTimeMillis();
//
//    static long wyhash64() {
//    	  wyhash64_x += 0x60bee2bee120fc15L;
//    	  long tmp =  wyhash64_x * 0xa3b195354a39b70dL;
//    	  final long  m1 = (tmp >> 64) ^ tmp;
//    	  tmp = m1 * 0x1b03738712fad5c9L;
//    	  final long m2 = (tmp >> 64) ^ tmp;
//    	  return m2;
//    }

	/**
	 * Draws a random number from a Poisson distribution using Knuth's algorithm.
	 */
	public static int drawFromPoisson(double lambda) {
	    double L = Math.exp(-lambda);
	    int k = 0;
	    double p = 1.0;
	    do {
	        k++;
	        p *= nextDouble();
	    } while (p > L);
	    return k - 1;
	}

	/**
	 * Draws a random state from a discrete categorical distribution.
	 */
	public static int drawFromCategorical(double[] probabilities) {
	    double u = nextDouble();
	    double cumulativeProbability = 0.0;
	    for (int i = 0; i < probabilities.length; i++) {
	        cumulativeProbability += probabilities[i];
	        if (u < cumulativeProbability) {
	            return i;
	        }
	    }
	    
	    Log.warning("Potential invalid probabilities (do not add up to 1)");
	    return probabilities.length - 1; // Should not be reached with valid probabilities
	}
}
