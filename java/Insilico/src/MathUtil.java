import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;


public class MathUtil {
	//  function to normalize vectors
	public static double[] norm(double[] x){
		double[] xnorm = new double[x.length];
		double sumx = 0;
		for(int i = 0; i < x.length; i++) sumx += x[i];
		if(sumx == 0){
			for(int i = 0 ; i < xnorm.length; i++) xnorm[i] = 1;
		}else{
			for(int i = 0; i < xnorm.length; i++) xnorm[i] = x[i] / sumx;
		}
		return(xnorm);
	}
	// function to sample inverse beta distribution.	
	public static double truncbeta(BetaDistribution beta, Random rand, 
			double min, double max){
		double value = min;
		double ymin = beta.cumulativeProbability(min);
		double ymax = beta.cumulativeProbability(max);
		// handling boundary case
		if(Math.abs(ymax - ymin) < 1e-8){
			//double mean = beta.getNumericalMean();
			return((max + min)/2.0);
			//return((mean < (max + min)/2.0 )? min : max);
		}
		value = beta.inverseCumulativeProbability(rand.nextDouble() * (ymax - ymin) + ymin);
		if(value == 0) 		System.out.printf("lower %.6f, upper %.6f, sampled 0\n", ymin, ymax);
		return(value);
	}
	// function to find min of selected elements in an array.
	public static double array_min(double[] array, ArrayList<Integer> location){
		double min = Double.MAX_VALUE;
		for(int i : location){
			if(array[i] < min) min = array[i];
		}
		return(min);
	}
	// function to find max of selected elements in an array.
	public static double array_max(double[] array,  ArrayList<Integer> location){
		double max = Double.MIN_VALUE;
		for(int i : location){
			if(array[i] > max) max = array[i];
		}
		return(max);
	}	
	
	// function to grab certain column of 2d array
	public static double[] grab2(double[][] matrix, int col){
		double[] out = new double[matrix.length];
		for(int i = 0; i < out.length; i++) out[i] = matrix[i][col];
		return(out);
	}
	// function to grab certain column of 2d array
	public static int[] grab2(int[][] matrix, int col){
			int[] out = new int[matrix.length];
			for(int i = 0; i < out.length; i++) out[i] = matrix[i][col];
			return(out);
	}

}
