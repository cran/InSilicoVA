import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;

import cern.jet.random.tdouble.Gamma;
import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;

//	import java.util.*;

public class InsilicoSampler {
	public int N;
	public int S;
	public int C;
	public int N_sub;
	public int N_level;
	public double[][] probbase; // S by C matrix
	public double[][] probbase_order;  // S by C matrix, value from 1 (I) to 15 (N)
	public double[] level_values; // vector of length N_level
	// key: 0 to C-1, second key: 1:N_level, value: which symptoms (0 to S-1)
	public HashMap<Integer, HashMap<Integer, ArrayList<Integer>>> probbase_level; 	
	public int[][] count_m; // S by C matrix, count of "Y"
	public int[][] count_m_all; // S by C matrix, count of "Y" and ""
	public int[] count_c; // C vector
	
	/*
	 * Initialization of Insilico Sampler
	 */
	
	public void initiate(int N, int S, int C, int N_sub, int N_level, int[] subpop,
			double[][] probbase, double[][] probbase_order, double[] level_values){
		this.N = N;
		this.S = S; 
		this.C = C;
		this.N_sub = N_sub;
		this.N_level = N_level;
		this.probbase = probbase;
		this.probbase_order = probbase_order;
		this.count_m = new int[S][C];
		this.count_m_all = new int[S][C];
		this.count_c = new int[C];
		this.level_values = level_values;
		this.probbase_level = new HashMap<Integer, HashMap<Integer, ArrayList<Integer>>>();
		this.levelize();
	}
	
	/*
	 * function to summarize the current levels of prob base.
	 */
	private void levelize(){
		for(int s = 0; s < S; s++){
			for(int c = 0; c < C; c++){
				int level = (int) this.probbase_order[s][c];
				if(this.probbase_level.get(c) == null) this.probbase_level.put(c,  new HashMap<Integer,ArrayList<Integer>>());
				if(this.probbase_level.get(c).get(level) == null) this.probbase_level.get(c).put(level, new ArrayList<Integer>());
				this.probbase_level.get(c).get(level).add(s);
			}
		}
	}

	/*
	 * function to count the current configurations for TruncBeta
	 * 
	 */
	public void CountCurrent(double[][] indic, int[] ynew){
		this.count_m = new int[S][C];
		this.count_m_all = new int[S][C];
		this.count_c = new int[C];

		for(int n = 0; n < this.N; n++){
			int c_current = ynew[n];
			this.count_c[c_current] += 1;
			for(int s = 0; s < this.S; s++){
				if(indic[n][s] == 1){
					this.count_m[s][c_current] += 1;
					this.count_m_all[s][c_current] += 1;
				}else if(indic[n][s] == 0){
					this.count_m_all[s][c_current] += 1;
				}
			}
		}
	}
	 
	/* 
	 * function to calculate p.nb with sub-population
	 * input: contains_missing, indic (N x S),csmf_sub (Nsub x C), subpop (N)
	 * output: nb (N by C)
	 */	
	public double[][] pnb(boolean contains_missing, 
			double[][] indic, double[][] csmf_sub, 
			int[] subpop){
		// initialize p.nb matrix to p.hat
		double[][] nb = new double[this.N][this.C];
		for(int n = 0; n < this.N; n++){
			for(int c = 0; c < this.C; c++){
				nb[n][c] =  csmf_sub[subpop[n]][c];
			}
		}
		// calculate posterior
		for(int n = 0 ; n < this.N; n++){
			ArrayList<Integer> nomissing = new ArrayList<Integer>();
			for(int s = 0 ; s < indic[n].length; s++){
				if(indic[n][s] >= 0) nomissing.add(s);
			}
			for(int c = 0; c < this.C; c++){
				for(int s : nomissing){
					nb[n][c] *= (indic[n][s]>0) ? (this.probbase[s][c]) : (1-this.probbase[s][c]);
				}
			}
			nb[n] = MathUtil.norm(nb[n]);
		}
		return(nb);
	}

	/*
	 * function to sample y, note the sampled y starts from 0!
	 */
	public int[] sampleY(double[][] pnb, Random rand){
		int[] y = new int[this.N];
		for(int n = 0; n < this.N; n++){
			double u = rand.nextDouble();
			double cum = 0;
			for(int c = 0; c< this.C; c++){
				cum += pnb[n][c];
				if(u < cum){
					y[n] = c;
					break;
				}
			}
		}
		return(y);
	}

	/*
	 * function to update theta
	 * input: jumprange, mu (C) sigma, theta (C), Y(C), N, jump.prop, rngN, rngU
	 */
	public double[] thetaBlockUpdate(double jumprange, double[] mu, double sigma2, double[] theta, 
			int[] Y, boolean jump_prop, Normal rngN, Random rand){
		double[] jump = new double[this.C];
		for(int c = 0; c < this.C; c++) jump[c] = jumprange;
		if(jump_prop){
			//todo: update proportional jump
		}
		double[] theta_new = new double[this.C];
		theta_new[0] = 1.0;
		double expsum = Math.exp(1.0);
		double expsum_new = Math.exp(1.0);
		for(int c = 1; c < this.C; c++){
			theta_new[c] = rngN.nextDouble(theta[c], jump[c]);
			expsum += Math.exp(theta[c]);
			expsum_new += Math.exp(theta_new[c]);
		}
		double logTrans = 0;
		for(int c = 0; c < this.C; c++){
			double diffquad = (theta_new[c] - mu[c])*(theta_new[c] - mu[c]) 
					- (theta[c]-mu[c])*(theta[c]-mu[c]);
			logTrans += Y[c] * (theta_new[c] - theta[c] - Math.log(expsum_new / expsum))
					 - 1/(2*sigma2) * diffquad;
		} 
		double u = Math.log(rand.nextDouble());
		if(logTrans >= u){
			//System.out.printf("+");
		}else{
			//System.out.printf("-");
			return(theta);
		}
		return(theta_new);
	}


	/*
	 * function to perform Truncated beta distribution for the whole conditional probability matrix.
	 * update cond.prob within the class. 
	 * Note: between causes comparison not performed here. Truncation only performed within same cause
	 * 
	 * input:
	 * rand: randome number generator
	 * prior_a: prior of alpha vector
	 * prior_b: prior of beta 
	 * trunc_max, trunc_Min: max/min of truncation
	 * used fields:
	 * probbase_order, level_values, probbase_level, count_m, count_m_all, count_c
	 *   key: 1 to C, second key: 1:N_level, value: which symptoms
	 *   HashMap<Integer, HashMap<Integer, ArrayList<Integer>>> probbase_level; 
	 */
	public void TruncBeta(Random rand, double[] prior_a, double prior_b, 
			double trunc_min, double trunc_max){
		double a = 0;
		double b = 0;
		double[][] new_probbase = new double[this.S][this.C];
		for(int s=0; s<this.S; s++){
			for(int c=0; c < this.C; c++){new_probbase[s][c] = this.probbase[s][c];}
		}
		// cause c
		for( int c = 0; c < C; c++){
			HashMap<Integer, ArrayList<Integer>> levels_under_c = probbase_level.get(c);
			double[] prob_under_c = MathUtil.grab2(this.probbase, c);
			double[] new_prob_under_c = new double[S];
			ArrayList<Integer> exist_levels_under_c = new ArrayList<Integer>();
			for(int l = 1; l <= this.N_level; l++){
				if(levels_under_c.get(l) != null) exist_levels_under_c.add(l);
			}
			// level l, ascending order
			for(int index = 0; index < exist_levels_under_c.size(); index++){
				int l_current = exist_levels_under_c.get(index);
				// symptom s
				for(int s : levels_under_c.get(l_current)){
					int count = count_m[s][c];
					int count_all = count_m_all[s][c];
					double lower = 0;
					double upper = 1;
					if (index == 0){
						int l_next = exist_levels_under_c.get(index + 1 );
						lower = MathUtil.array_max(prob_under_c, levels_under_c.get(l_next));
						lower = Math.max(lower, trunc_min);
						upper = trunc_max;
					}else if(index == exist_levels_under_c.size() - 1){
						lower = trunc_min;
						int l_prev = exist_levels_under_c.get(index-1);
						upper = MathUtil.array_min(new_prob_under_c , levels_under_c.get(l_prev));
						upper = Math.min(upper, trunc_max);
					}else{
						int l_next = exist_levels_under_c.get(index + 1 );
						lower = MathUtil.array_max(prob_under_c, levels_under_c.get(l_next));
						lower = Math.max(lower, trunc_min);
						int l_prev = exist_levels_under_c.get(index - 1);
						upper = MathUtil.array_min(new_prob_under_c , levels_under_c.get(l_prev));
						upper = Math.min(upper,  trunc_max);
					}
					if(lower >= upper){new_prob_under_c[s] = upper;
					}else{
						a = prior_a[l_current-1] + count;
						b = prior_b + count_all - a;
						BetaDistribution beta = new BetaDistribution(a, b,  1e-10);
						new_prob_under_c[s] = MathUtil.truncbeta(beta, rand, lower, upper);
					}
				}
			}
			// update this column of probbase
			for(int s = 0; s < this.S; s++){this.probbase[s][c] = new_prob_under_c[s];}
		}

	}

	/*
	 * function to perform Truncated beta distribution for the probbase table.
	 * update cond.prob within the class. Return table vector.
	 * 
	 * input:
	 * rand: randome number generator
	 * prior_a: prior of alpha vector
	 * prior_b: prior of beta 
	 * last: current probbase
	 * trunc_max, trunc_Min: max/min of truncation
	 * used fields:
	 * probbase_order, level_values, probbase_level, count_m, count_m_all, count_c
	 *   key: 1 to C, second key: 1:N_level, value: which symptoms
	 *   HashMap<Integer, HashMap<Integer, ArrayList<Integer>>> probbase_level; 
	 */
	public void TruncBeta_pool(Random rand, double[] prior_a, double prior_b, 
			double trunc_min, double trunc_max){
		double a = 0;
		double b = 0;
		double[][] new_probbase = new double[this.S][this.C];
		for(int s=0; s<this.S; s++){
			for(int c=0; c < this.C; c++){new_probbase[s][c] = this.probbase[s][c];}
		}
		// assume all levels exist in data, if not, modify the R codes before calling
		double[] new_level_values = new double[this.N_level];
		for(int l = 1; l <= this.N_level; l++){
			int count = 0; 
			int count_all = 0;
			// count appearances
			for(int c = 0; c<this.C; c++){
				if(this.probbase_level.get(c).get(l) != null){
					for(int s : this.probbase_level.get(c).get(l)){
						count += count_m[s][c];
						count_all += count_m_all[s][c];
					}
				}
			} 
			double lower = 0;
			double upper = 1;
			// note l starts from 1, need to extract an extra 1
			if(l == 1){
				lower = Math.max(this.level_values[l], trunc_min);
				upper = trunc_max;
			}else if(l == this.N_level){
				lower = trunc_min;
				upper = Math.min(new_level_values[l-2], trunc_max);
				//upper = Math.min(this.level_values[l-2], upper);
			}else{
				lower = Math.max(this.level_values[l], trunc_min);
				upper =  Math.min(new_level_values[l-2], trunc_max);
				//upper = Math.min(this.level_values[l-2], upper);
			}
			if(lower >= upper){new_level_values[l-1] = upper;
			}else{
				a = prior_a[l-1] + count;
				b = prior_b + count_all - a;
				BetaDistribution beta = new BetaDistribution(a, b,  1e-10);
				new_level_values[l-1] = MathUtil.truncbeta(beta, rand, lower, upper);
			}	
		}
		this.level_values = new_level_values;
		for(int s = 0; s < this.S; s++){
			for(int c = 0; c < this.C; c++){
				this.probbase[s][c] = this.level_values[(int) this.probbase_order[s][c] - 1];
			}
		}
		return;
	}

	/*
	 * Main function to perform insilico sampling
	 * 
	 * key parameters:
	 * indic: N by S matrix
	 * subpop: length N vector of sub-population assignment
	 *         Note: should start from 0. All 0 if no sub-population exist.
	 * contains_missing: if there are missing
	 * pool: if estimating the probbase table only
	 * seed, N_gibbs, thin: integers
	 * mu: vector initialized in R
	 * sigma2: value initialized in R
	 */
	public static double[][] Fit(int N, int S, int C, int N_sub, int N_level, 
			double[][] probbase, double[][] probbase_order, double[] level_values, 
			double[] prior_a, double prior_b, double jumprange, double trunc_min, double trunc_max,
			double[][] indic, int[] subpop, int contains_missing, int pool,
			int seed, int N_gibbs, int burn, int thin, 
			double[] mu, double sigma2, boolean this_is_Unix, boolean useProbbase){
//			public static void main(String[] args){
//				int N = 5; 
//				int S = 3; 
//				int C = 2; 
//				int N_sub = 2; 
//				int N_level = 4; 
//				double[][] probbase = {{0.8, 0.1}, {0.2, 0.7}, {0.7,0.7}}; 
//				double[][] probbase_order = {{1, 4}, {3, 2}, {2,2}}; 
//				double[] level_values = {0.9, 0.8, 0.7,0.2}; 
//				double[] prior_a = {0.9, 0.8, 0.7,0.2}; 
//				double prior_b = 20; 
//				double jumprange = 10; 
//				double trunc_min = 0.01; 
//				double trunc_max = 0.99;
//				double[][] indic = {{1,0,-1},{0,0,1},{1,1,0},{1,0,1},{0,1,1}}; 
//				int[] subpop = {0,0,0,1,1}; 
//				int contains_missing = 1; 
//				int pool = 1;
//				int seed = 1; 
//				int burn = 0;
//				int N_gibbs = 500; 
//				int thin = 1; 
//				double[] mu = {0.2,0.2}; 
//				double sigma2 = 1;
//		        boolean this_is_Unix = true;
//				boolean useProbbase = false;
		
		
		InsilicoSampler insilico = new InsilicoSampler();
		insilico.initiate(N, S, C, N_sub, N_level, subpop, probbase, probbase_order, level_values);
		System.out.printf("Insilico Sampler initiated, %d iterations to sample\n", N_gibbs);
		DoubleRandomEngine rngEngine=new DoubleMersenneTwister(seed);
		Normal rngN=new Normal(0.0,1.0,rngEngine);
		Gamma rngG=new Gamma(1.0,1.0,rngEngine);
		Random rand = new Random(seed);

		// allocate values
		int N_thin = ((int) ((N_gibbs-burn) / (thin + 0.0)));
		int n_report = Math.max(N_gibbs / 20, 100);
		if(N_gibbs < 200) n_report = 50;
		
		//double[][][] mu_gibbs = new double[N_thin][N_sub][C];
		//double[][] sigma2_gibbs = new double[N_thin][N_sub];
		double[][][] probbase_gibbs = new double[N_thin][S][C];
		double[][] levels_gibbs = new double[N_thin][N_level];
		//double[][][] theta_gibbs = new double[N_thin][N_sub][C];
		double[][][] p_gibbs = new double[N_thin][N_sub][C];
		double[][][] pnb_gibbs = new double[N_thin][N][C];

		int[] naccept = new int[N_sub];
		double[][] mu_now = new double[N_sub][C];
		double[] sigma2_now = new double[N_sub];
		double[][] theta_now = new double[N_sub][C];
		double[][] p_now = new double[N_sub][C];

		// initialize values
		for(int sub  = 0; sub < N_sub; sub++){
			mu_now[sub] = mu;
			sigma2_now[sub] = sigma2;
			theta_now[sub][0] = 1;
			double expsum = Math.exp(1.0);
			for(int c = 1; c < C; c++){
				theta_now[sub][c] = Math.log(rand.nextDouble() * 100.0);
				expsum += Math.exp(theta_now[sub][c]);
			}
			for(int c = 0; c < C; c++){
				p_now[sub][c] = Math.exp(theta_now[sub][c]) / expsum;
			}
		}

		// first time pnb calculation, note it is inverse of R version
		double[][] pnb = new double[N][C];
		pnb = insilico.pnb((contains_missing == 1), indic, p_now, subpop);


		// start loop
		long start = System.currentTimeMillis();
	   
		ProgressPopup popup = new ProgressPopup(this_is_Unix, N_gibbs);
		
		for(int k = 0; k < N_gibbs; k++){
			if(!this_is_Unix){
				popup.update(k);
			} 
			// sample new y vector
			int[] y_new = insilico.sampleY(pnb, rand);
			// count the appearance of each cause
			int[][] Y = new int[N_sub][C];
			for(int n = 0; n < N; n++){Y[subpop[n]][y_new[n]] += 1;}

			for(int sub = 0; sub < N_sub; sub++){
				// sample mu
				double mu_mean = 0;
				for(int c = 0; c<C; c++) mu_mean += theta_now[sub][c];
				mu_mean = mu_mean / (C+0.0);
				mu_mean = rngN.nextDouble(mu_mean, Math.sqrt(sigma2_now[sub]/ (C+0.0)));
				for(int c = 0; c < C; c++){
					mu_now[sub][c] = mu_mean;
				}

				// sample sigma2
				double shape = (C-1.0)/2;
				double rate2 = 0;
				for(int c = 0 ; c < C; c++) rate2 += Math.pow(theta_now[sub][c] - mu_now[sub][c], 2);
				sigma2_now[sub] = 1 / rngG.nextDouble( shape, rate2/2 );
				// sample theta
				double[] theta_prev = theta_now[sub];
				theta_now[sub] = insilico.thetaBlockUpdate(jumprange, mu_now[sub], 
						 sigma2_now[sub], theta_prev, Y[sub], false, rngN, rand);	
				if(theta_now[sub][1] != theta_prev[1]) naccept[sub] += 1; 
				
				// calculate phat
				double expsum = 0.0;
				for(int c = 0; c < C; c++){
					expsum += Math.exp(theta_now[sub][c]);
				}
				for(int c = 0; c < C; c++){
					p_now[sub][c] = Math.exp(theta_now[sub][c]) / expsum;
				}
			}
			
 			if(useProbbase){
 				// skip the update step for prob base
 			}else{
 				insilico.CountCurrent(indic,  y_new);
 				if(pool == 1){
 					insilico.TruncBeta_pool(rand, prior_a, prior_b, trunc_min, trunc_max);				
 				}else{
 					insilico.TruncBeta(rand, prior_a, prior_b, trunc_min, trunc_max);	
 				}
 			}
			
			pnb = insilico.pnb((contains_missing==1), indic, p_now, subpop);

			// format output message
			if(k % 10 == 0) System.out.printf(".");
			if(k % n_report == 0 & k != 0){
				// output for Unix system
				long now   = System.currentTimeMillis();
				String message = String.format("\nIteration: %d \n", k);
				for(int sub = 0; sub < N_sub; sub++){
					double ratio = naccept[sub] / (k+0.0);
					message += String.format("Sub-population %d acceptance ratio: %.2f \n", sub, ratio);					
				}
				System.out.printf(message);
				System.out.printf("%.2fmin elapsed, expected %.2fmin remaining \n", 
						(double) (now - start)/1000.0/60.0, 
						(double) (now - start)/1000.0/60.0/(k+0.0) * (N_gibbs - k));
				
				// output for windows pop up window
				if(!this_is_Unix){
					popup.message(k, message);
				}
			}
						

			// note this condition includes the first iteration
			if(k >= burn && (k-burn+1) % thin == 0){
				int save = ((int) ((k-burn+1) / thin)) - 1;
				//mu_gibbs[save] = mu_now;
				//sigma2_gibbs[save] = sigma2_now;
				//theta_gibbs[save] = theta_now;
				for(int d1 = 0; d1 < N; d1++){
					for(int d2 = 0; d2 < C; d2++){
						pnb_gibbs[save][d1][d2] = pnb[d1][d2];
					}
				}
				for(int d1 = 0; d1 < N_sub; d1++){
					for(int d2 = 0; d2 < C; d2++){
						p_gibbs[save][d1][d2] = p_now[d1][d2];		
					}
				}
				if(pool == 1){
					for(int d1 = 0; d1 < N_level; d1++){
						levels_gibbs[save][d1] =  insilico.level_values[d1];
					}
				}else{	
					for(int d1 = 0; d1 < S; d1++){
						for(int d2 = 0; d2 < C; d2++)
							probbase_gibbs[save][d1][d2] = insilico.probbase[d1][d2];
					}
				}
			}
		}
		// close windows pop up window
		if(!this_is_Unix) { popup.close();}	
		/*
		 * format output from here
		 */
		System.out.println("\nOverall acceptance ratio:");
		for(int sub = 0; sub < N_sub; sub++){
			double ratio = naccept[sub] / (N_gibbs+0.0);
			System.out.printf("Sub-population %d : %.4f \n", sub, ratio);
		}	
		int N_out = 0;
		if(pool == 1){
			N_out =  N_sub * C + N * C + N_level;
		}else{
			N_out =  N_sub * C + N * C + S * C;
		}

		/*      
		 * 	 Parameters for output:		
		 *		double[][][] p_gibbs = new double[N_thin][N_sub][C];
		 *		int[][] y_gibbs = new int[N_thin][N];
		 *		double[][][] probbase_gibbs = new double[N_thin][S][C];
		 *		double[][] levels_gibbs = new double[N_thin][N_level];
		 */
		double[][] parameters = new double[N_thin][N_out];
		for(int k = 0; k < N_thin; k++){
			int counter = 0;
			// save p_gibbs
			for(int sub = 0; sub < N_sub; sub++){
				for(int c = 0; c < C; c++){
					parameters[k][counter] = p_gibbs[k][sub][c];
					counter += 1;
				}
			}
			// save pnb_gibbs
			for(int n = 0; n < N; n++){
				for(int c = 0; c < C; c++){
					parameters[k][counter] = pnb_gibbs[k][n][c]; 
					counter += 1;					
				}
			}
			if(pool != 1){
				for(int s = 0; s < S; s++){
					for(int c = 0; c < C; c++){
						parameters[k][counter] = probbase_gibbs[k][s][c];
						counter += 1;
					}
				}
			}else{
				for(int l = 1; l <= N_level; l++){parameters[k][counter] = levels_gibbs[k][l-1]; counter+=1;}
			}
		}
		System.out.println("Organizing output, might take a moment...");
		
//		return;
		return(parameters);
	}

}