/**
 * Created by zehangli on 9/19/15.
 */
package sampler;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import cern.colt.Arrays;
import utils.MathUtil;
import utils.ProgressPopup;
import utils.ProgressBar;

import org.apache.commons.math3.distribution.BetaDistribution;

import cern.jet.random.tdouble.Gamma;
import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;


public class InsilicoSampler2 {
    // number of death
    public int N;
    // number of symptoms
    public int S;
    // number of causes
    public int C;
    // number of sub populations
    public int N_sub;
    // number of levels in the probbase table
    public int N_level;
    // S by C matrix of probbase
    public double[][] probbase;
    // S by C matrix of probbase in order index, value from 1 (I) to 15 (N)
    public double[][] probbase_order;
    // vector of length N_level, prior of the level values
    public double[] level_values;
    // HashMap that tells which symptoms are at which level in each causes
    // Top layer:     key = 0 to C-1
    // Second layer:  key = 1 to N_level
    // Value Set:     which symptoms (0 to S-1)
    public HashMap<Integer, HashMap<Integer, ArrayList<Integer>>> probbase_level;
    // S by C matrix, counts of "Y" (yes)
    public int[][] count_m;
    // S by C matrix, counts of "Y" and "" (yes and no, not missing)
    public int[][] count_m_all;
    // length C vector, counts of each causes
    public int[] count_c;

    // stuff related to physician coding
    // number of physician coded categories
    public int C_phy;
    // length C conversion vector
    public double[] broader;
    // HashMap that tells which causes go into each category
    // Key: category, Value: set of causes
    public HashMap<Integer, HashSet<Integer>> catmap;

    /*
     * Initialization of InSilico Sampler
     */
    public void initiate(int N, int S, int C, int N_sub, int N_level, int[] subpop,
                         double[][] probbase, double[][] probbase_order, double[] level_values, int pool) {
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
        this.levelize(pool);
    }
    /*
     * Initialization with physician coding
     */
    public void initiate(int N, int S, int C, int N_sub, int N_level, int[] subpop,
                         double[][] probbase, double[][] probbase_order, double[] level_values, int pool,
                         int C_phy, double[] broader){
        this.initiate(N, S, C, N_sub, N_level, subpop, probbase, probbase_order, level_values, pool);
        // added variables
        this.C_phy = C_phy;
        this.broader = broader;
        this.initiateCatMap();
    }

    /*
     * Initialize Category Mapping
     */
    public void initiateCatMap(){
        this.catmap = new HashMap<Integer, HashSet<Integer>>();
        // unknown causes
        this.catmap.put(0, new HashSet<Integer>());
        for(int j = 1; j < this.broader.length; j++) {
            this.catmap.get(0).add(j);
        }
        // actual causes
        for(int i = 1; i < this.C_phy; i++){
            this.catmap.put(i, new HashSet<Integer>());
            for(int j = 1; j < this.broader.length; j++){
                if( (int) this.broader[j] == i){
                    this.catmap.get(i).add(j);
                }
            }
        }
    }
    /*
     * function to summarize the current levels of prob base.
     */
    private void levelize(int pool){
        // if pool = 0, doesn't matter which
        // if pool = 1, count level by cause
        if(pool == 0 | pool == 1){
            // loop over all s and c combination
            for(int s = 0; s < S; s++) {
                for (int c = 0; c < C; c++) {
                    // get the level of this s-c combination
                    int level = (int) this.probbase_order[s][c];
                    // initialize if this cause has not been initialized in probbase_level yet
                    if (this.probbase_level.get(c) == null)
                        this.probbase_level.put(c, new HashMap<Integer, ArrayList<Integer>>());
                    // initialize if this level under this cause has not been initialized
                    if (this.probbase_level.get(c).get(level) == null)
                        this.probbase_level.get(c).put(level, new ArrayList<Integer>());
                    // save the cause-level-symptom combination
                    this.probbase_level.get(c).get(level).add(s);
                }
            }
        // if pool = 2, count level by symptom
        }else if(pool == 2){
            // loop over all s and c combination
            for(int s = 0; s < S; s++) {
                for (int c = 0; c < C; c++) {
                    // get the level of this s-c combination
                    int level = (int) this.probbase_order[s][c];
                    // initialize if this cause has not been initialized in probbase_level yet
                    if (this.probbase_level.get(s) == null)
                        this.probbase_level.put(s, new HashMap<Integer, ArrayList<Integer>>());
                    // initialize if this level under this symptom has not been initialized
                    if (this.probbase_level.get(s).get(level) == null)
                        this.probbase_level.get(s).put(level, new ArrayList<Integer>());
                    // save the cause-level-symptom combination
                    this.probbase_level.get(s).get(level).add(c);
                }
            }
        }
    }

    /*
     * function to count the current configurations for TruncBeta
     *
     */
    public void CountCurrent(double[][] indic, int[] ynew){
        // counts of yes
        this.count_m = new int[S][C];
        // counts of yes and no
        this.count_m_all = new int[S][C];
        // counts of causes
        this.count_c = new int[C];

        // loop over all individuals
        for(int n = 0; n < this.N; n++){
            // cause of this death
            int c_current = ynew[n];
            // add to total counts
            this.count_c[c_current] += 1;
            // loop over all symptoms
            for(int s = 0; s < this.S; s++){
                // add to counts of this symptom's appearance
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
     * some key dimensions: indic (N x S),csmf_sub (Nsub x C), subpop (N)
     * output: nb (N by C)
     */
    public double[][] pnb(boolean contains_missing,
                          double[][] indic, double[][] csmf_sub,
                          int[] subpop, int[][] zero_matrix){
        // initialize p.nb matrix to p.hat
        double[][] nb = new double[this.N][this.C];
        // pnb_{nc} <- csmf_{subpop_n, c}
        for(int n = 0; n < this.N; n++){
            for(int c = 0; c < this.C; c++){
                nb[n][c] =  csmf_sub[subpop[n]][c] * zero_matrix[n][c];
            }
        }
        // calculate posterior
        for(int n = 0 ; n < this.N; n++){
            // find which symptoms are not missing for this death
            ArrayList<Integer> nomissing = new ArrayList<Integer>();
            for(int s = 0 ; s < indic[n].length; s++){
                if(indic[n][s] >= 0) nomissing.add(s);
            }
            // loop over cause-symptoms combination to calculate naive bayes prob
            for(int c = 0; c < this.C; c++){
                for(int s : nomissing){
                    nb[n][c] *= (indic[n][s]>0) ? (this.probbase[s][c]) : (1-this.probbase[s][c]);
                }
            }
            // normalization
            nb[n] = MathUtil.norm(nb[n]);
        }
        return(nb);
    }

    /*
     * function to calculate p.nb with sub-population and physician coding
     * some key dimensions: indic (N x S),csmf_sub (Nsub x C), subpop (N), assignment(N)
     * output: nb (N by C)
     */
    public double[][] pnb(boolean contains_missing,
                          double[][] indic, double[][] csmf_sub,
                          int[] subpop, int[] g_new, int[][] zero_matrix){
        // initialize p.nb matrix to p.hat
        double[][] nb = new double[this.N][this.C];
        // pnb_{nc} <- csmf_{subpop_n, c}

        // sample a cause group first, using the same as sample new Y
        //default to zero, update those with physician codes, normalize at the end
        for(int n = 0; n < this.N; n++){
            for(int c : this.catmap.get(g_new[n])){
                nb[n][c] =  csmf_sub[subpop[n]][c] * zero_matrix[n][c];
            }
        }

        // calculate posterior
        for(int n = 0 ; n < this.N; n++){
            // find which symptoms are not missing for this death
            ArrayList<Integer> nomissing = new ArrayList<Integer>();
            for(int s = 0 ; s < indic[n].length; s++){
                if(indic[n][s] >= 0) nomissing.add(s);
            }
            // loop over cause-symptoms combination to calculate naive bayes prob
            for(int c : this.catmap.get(g_new[n])){
                for(int s : nomissing){
                    nb[n][c] *= (indic[n][s]>0) ? (this.probbase[s][c]) : (1-this.probbase[s][c]);
                }
            }
            // normalization
            nb[n] = MathUtil.norm(nb[n]);
        }
        return(nb);
    }

    /*
     * function to sample a multinomial, note the sampled y starts from 0!
     */
    public int[] sampleY(double[][] pnb, Random rand){
        int[] y = new int[pnb.length];
        // loop over every death
        for(int n = 0; n < pnb.length; n++){
            // naive implementation of categorical distribution
            double u = rand.nextDouble();
            double cum = 0;
            for(int c = 0; c< pnb[n].length; c++){
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
     * input (dimensions): jumprange, mu (C) sigma, theta (C), Y(C), N, jump.prop(NOT USING), rngN, rngU
     * new input: zero_vector: define which causes are set to zero
     */
    public double[] thetaBlockUpdate(double jumprange, double[] mu, double sigma2, double[] theta,
                                     int[] Y, boolean jump_prop, Normal rngN, Random rand, int[] zero_vector){
        // initialize jump range, use the same value for all causes
        double[] jump = new double[this.C];
        for(int c = 0; c < this.C; c++) jump[c] = jumprange;

        if(jump_prop){
            // not using
        }

        // initialize new theta
        double[] theta_new = new double[this.C];
        // fix first theta (not removed)
        int fix = 0;
        for(int i = 0; i < zero_vector.length; i++){
            if(zero_vector[i] != 0){
                fix = i;
                break;
            }
        }
        theta_new[fix] = 1.0;
        // calculate sum(exp(theta)), the normalizing constant
        double expsum = Math.exp(1.0);
        double expsum_new = Math.exp(1.0);

        // sample new theta proposal
        for(int c = fix + 1; c < this.C; c++){
            theta_new[c] = rngN.nextDouble(theta[c], jump[c]) * zero_vector[c];
            expsum += Math.exp(theta[c]);
            expsum_new += Math.exp(theta_new[c]);
        }

        // calculate likelihood
        double logTrans = 0;
        // l = Y * (theta_new - theta - log(expsum_new/expsum)) - ((theta_new - mu)^2 - (theta - mu)^2 ) / (2*sigma^2)
        for(int c = 0; c < this.C; c++){
            if(zero_vector[c] > 0){
                double diffquad = (theta_new[c] - mu[c])*(theta_new[c] - mu[c])
                        - (theta[c]-mu[c])*(theta[c]-mu[c]);
                logTrans += Y[c] * (theta_new[c] - theta[c] - Math.log(expsum_new / expsum))
                        - 1/(2*sigma2) * diffquad;
            }
        }

        // accept or reject
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
    * Update Sep 30, 2015: truncated beta for new SCI
    * function to perform Truncated beta distribution for the whole conditional probability matrix.
    * update cond.prob within the class.
    * Note: between symptom comparison not performed here. Truncation only performed within same symptom
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
    public void TruncBeta2(Random rand, double[] prior_a, double prior_b,
                          double trunc_min, double trunc_max){
        double a = 0;
        double b = 0;
        // create a new transpose probbase matrix
        double[][] new_probbase = new double[this.C][this.S];
        for(int s=0; s<this.S; s++){
            for(int c=0; c < this.C; c++){new_probbase[c][s] = this.probbase[s][c];}
        }
        // loop over causes c
        for( int s = 0; s < this.S; s++){
            // find which level-symptom combinations under this cause
            HashMap<Integer, ArrayList<Integer>> levels_under_s = this.probbase_level.get(s);
            // get the conditional probabilities of each symptoms given cause c
            double[] prob_under_s = this.probbase[s];
            // create new instance of this vector
            double[] new_prob_under_s = new double[this.C];
            // find the list of all levels present under cause c
            ArrayList<Integer> exist_levels_under_s = new ArrayList<Integer>();
            for(int l = 1; l <= this.N_level; l++){
                if(levels_under_s.get(l) != null) exist_levels_under_s.add(l);
            }
            // loop over level l, in ascending order
            for(int index = 0; index < exist_levels_under_s.size(); index++){
                int l_current = exist_levels_under_s.get(index);
                // loop over symptoms s in the level l
                for(int c : levels_under_s.get(l_current)){
                    // find the counts of this symptom
                    int count = count_m[s][c];
                    int count_all = count_m_all[s][c];
                    double lower = 0;
                    double upper = 1;
                    // if the highest level
                    if (index == 0){
                        // find which level is next
                        int l_next = exist_levels_under_s.get(index + 1 );
                        // find the max of those symptoms
                        lower = MathUtil.array_max(prob_under_s, levels_under_s.get(l_next));
                        // make sure not lower than lower bound
                        lower = Math.max(lower, trunc_min);
                        upper = trunc_max;
                        // if the lowest level
                    }else if(index == exist_levels_under_s.size() - 1){
                        lower = trunc_min;
                        int l_prev = exist_levels_under_s.get(index-1);
                        upper = MathUtil.array_min(new_prob_under_s, levels_under_s.get(l_prev));
                        upper = Math.min(upper, trunc_max);
                        // if in the middle
                    }else{
                        int l_next = exist_levels_under_s.get(index + 1 );
                        lower = MathUtil.array_max(prob_under_s, levels_under_s.get(l_next));
                        lower = Math.max(lower, trunc_min);
                        int l_prev = exist_levels_under_s.get(index - 1);
                        upper = MathUtil.array_min(new_prob_under_s, levels_under_s.get(l_prev));
                        upper = Math.min(upper,  trunc_max);
                    }
                    // if range is invalid, use higher case
                    if(lower >= upper){
                        new_prob_under_s[c] = upper;
                    }else{
                        // find the beta distribution parameters, note level starting from 1
                        a = prior_a[l_current-1] + count;
                        b = prior_b + count_all - a;
                        BetaDistribution beta = new BetaDistribution(a, b,  1e-10);
                        new_prob_under_s[c] = MathUtil.truncbeta(beta, rand, lower, upper);
                    }
                }
            }
            // update this column of probbase
            for(int c = 0; c < this.C; c++){this.probbase[s][c] = new_prob_under_s[c];}
        }
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
        // create a new probbase matrix
        double[][] new_probbase = new double[this.S][this.C];
        for(int s=0; s<this.S; s++){
            for(int c=0; c < this.C; c++){new_probbase[s][c] = this.probbase[s][c];}
        }
        // loop over causes c
        for( int c = 0; c < C; c++){
            // find which level-symptom combinations under this cause
            HashMap<Integer, ArrayList<Integer>> levels_under_c = this.probbase_level.get(c);
            // get the conditional probabilities of each symptoms given cause c
            double[] prob_under_c = MathUtil.grab2(this.probbase, c);
            // create new instance of this vector
            double[] new_prob_under_c = new double[this.S];
            // find the list of all levels present under cause c
            ArrayList<Integer> exist_levels_under_c = new ArrayList<Integer>();
            for(int l = 1; l <= this.N_level; l++){
                if(levels_under_c.get(l) != null) exist_levels_under_c.add(l);
            }
            // loop over level l, in ascending order
            for(int index = 0; index < exist_levels_under_c.size(); index++){
                int l_current = exist_levels_under_c.get(index);
                // loop over symptoms s in the level l
                for(int s : levels_under_c.get(l_current)){
                    // find the counts of this symptom
                    int count = count_m[s][c];
                    int count_all = count_m_all[s][c];
                    double lower = 0;
                    double upper = 1;
                    // if the highest level
                    if (index == 0){
                        // find which level is next
                        int l_next = exist_levels_under_c.get(index + 1 );
                        // find the max of those symptoms
                        lower = MathUtil.array_max(prob_under_c, levels_under_c.get(l_next));
                        // make sure not lower than lower bound
                        lower = Math.max(lower, trunc_min);
                        upper = trunc_max;
                        // if the lowest level
                    }else if(index == exist_levels_under_c.size() - 1){
                        lower = trunc_min;
                        int l_prev = exist_levels_under_c.get(index-1);
                        upper = MathUtil.array_min(new_prob_under_c, levels_under_c.get(l_prev));
                        upper = Math.min(upper, trunc_max);
                        // if in the middle
                    }else{
                        int l_next = exist_levels_under_c.get(index + 1 );
                        lower = MathUtil.array_max(prob_under_c, levels_under_c.get(l_next));
                        lower = Math.max(lower, trunc_min);
                        int l_prev = exist_levels_under_c.get(index - 1);
                        upper = MathUtil.array_min(new_prob_under_c, levels_under_c.get(l_prev));
                        upper = Math.min(upper,  trunc_max);
                    }
                    // if range is invalid, use higher case
                    if(lower >= upper){
                        new_prob_under_c[s] = upper;
                    }else{
                        // find the beta distribution parameters, note level starting from 1
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
        // initialize a new probbase copy
        double[][] new_probbase = new double[this.S][this.C];
        for(int s=0; s<this.S; s++){
            for(int c=0; c < this.C; c++){new_probbase[s][c] = this.probbase[s][c];}
        }
        // assume all levels exist in data, if not, need to modify the R codes before calling
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
            // note l starts from 1, level_values have index starting from 0
            if(l == 1){
                // lower bound is the max of this next level
                lower = Math.max(this.level_values[l], trunc_min);
                upper = trunc_max;
            }else if(l == this.N_level){
                lower = trunc_min;
                // upper bound is the min of previous level
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
     *       update Sep 30, 2015: pool = 0, only table; pool = 1, by cause; pool = 2, compare only within each symptom
     * seed, N_gibbs, thin: integers
     * mu: vector initialized in R
     * sigma2: value initialized in R
     */
    public static double[] Fit(int[] dimensions,
                               double[][] probbase, double[][] probbase_order, double[] level_values,
                               double[] prior_a, double prior_b, double jumprange, double trunc_min, double trunc_max,
                               double[][] indic, int[] subpop, int contains_missing, int pool,
                               int seed, int N_gibbs, int burn, int thin,
                               double[] mu, double sigma2, boolean this_is_Unix, boolean useProbbase,
                               boolean isAdded,
                               double[][] mu_continue, double[] sigma2_continue, double[][] theta_continue,
                               int C_phy, double[] broader, double[][] assignment, int[][] impossible){
//			public static void main(String[] args){
//				int N = 5;
//				int S = 3;
//				int C = 2;
//				int N_sub = 2;
//				int N_level = 4;
//				double[][] probbase = {{0.8, 0.1}, {0.2, 0.7}, {0.7,0.7}};
//				double[][] probbase_order = {{1, 4}, {3, 2}, {2,3}};
//				double[] level_values = {0.9, 0.8, 0.7,0.2};
//				double[] prior_a = {0.9, 0.8, 0.7,0.2};
//				double prior_b = 20;
//				double jumprange = 10;
//				double trunc_min = 0.01;
//				double trunc_max = 0.99;
//				double[][] indic = {{1,0,-1},{0,0,1},{1,1,0},{1,0,1},{0,1,1}};
//				int[] subpop = {0,0,0,1,1};
//				int contains_missing = 1;
//				int pool = 2;
//				int seed = 1;
//				int burn = 0;
//				int N_gibbs = 500;
//				int thin = 1;
//				double[] mu = {0.2,0.2};
//				double sigma2 = 1;
//		        boolean this_is_Unix = true;
//				boolean useProbbase = false;
//                boolean isAdded = false;
//                double[][] mu_continue = {{0,0}, {0,0}};
//                double[] sigma2_continue = {0,0};
//                double[][] theta_continue = {{0,0}, {0,0}};
//                int C_phy = 2;
//                double[] broader = {1,1};
//                double[][] assignment = {{1,0}, {0,1}, {0.5,0.5}, {1,0}, {0,1}};

        // initialization
        int N = dimensions[0];
        int S = dimensions[1];
        int C = dimensions[2];
        int N_sub = dimensions[3];
        int N_level = dimensions[4];

        InsilicoSampler2 insilico = new InsilicoSampler2();
        boolean withPhy = (C_phy > 1);

        if(!withPhy){
            insilico.initiate(N, S, C, N_sub, N_level, subpop, probbase, probbase_order, level_values, pool);
        }else{
            insilico.initiate(N, S, C, N_sub, N_level, subpop, probbase, probbase_order, level_values, pool, C_phy,
                    broader);
        }

        System.out.printf("InSilico Sampler initiated, %d iterations to sample\n", N_gibbs);
        // list of random number generators to use
        DoubleRandomEngine rngEngine=new DoubleMersenneTwister(seed);
        Normal rngN=new Normal(0.0,1.0,rngEngine);
        Gamma rngG=new Gamma(1.0,1.0,rngEngine);
        Random rand = new Random(seed);

        // calculate the dimension of values to save, and the interval of report (print message on screen)
        int N_thin = ((int) ((N_gibbs-burn) / (thin + 0.0)));
        int n_report = Math.max(N_gibbs / 20, 100);
        if(N_gibbs < 200) n_report = 50;

        // probbase at each thinned iteration
        double[][][] probbase_gibbs = new double[N_thin][S][C];
        // probbase levels at each thinned iteration
        double[][] levels_gibbs = new double[N_thin][N_level];
        // csmf at each thinned iteration
        double[][][] p_gibbs = new double[N_thin][N_sub][C];
        // individual probability at each thinned iteration
        double[][] pnb_mean = new double[N][C];
        // number of acceptance
        int[] naccept = new int[N_sub];
        // current mu
        double[][] mu_now = new double[N_sub][C];
        // current sigma^2
        double[] sigma2_now = new double[N_sub];
        // current theta
        double[][] theta_now = new double[N_sub][C];
        // current csmf
        double[][] p_now = new double[N_sub][C];

        // if the chain is the original one, initialize randomly
        if(!isAdded){
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
        }else{
            // if the chain is to continue from a previous one, initialize with the _now values imported
            for(int sub  = 0; sub < N_sub; sub++){
                mu_now[sub] = mu_continue[sub];
                sigma2_now[sub] = sigma2_continue[sub];
                theta_now[sub] = theta_continue[sub];
                // recalculate p from theta
                double expsum = Math.exp(1.0);
                for(int c = 1; c < C; c++){
                    expsum += Math.exp(theta_now[sub][c]);
                }
                for(int c = 0; c < C; c++){
                    p_now[sub][c] = Math.exp(theta_now[sub][c]) / expsum;
                }
            }

        }

        // check impossible causes?
        boolean check_impossible = impossible[0].length == 2;
        int[][] zero_matrix = new int[N][C];
        for(int i = 0; i < N; i++){
            for(int j = 0; j < C; j++){
                zero_matrix[i][j] = 1;
            }
        }
        if(check_impossible){
            for(int i = 0; i < N; i++){
                for(int k = 0; k < impossible.length; k++){
                    if(indic[i][impossible[k][1] - 1] == 1){
                        zero_matrix[i][impossible[k][0] - 1] = 0;
                    }
                }
            }
        }
        // check if specific causes are impossible for a whole subpopulation
        int[][] zero_group_matrix = new int[N_sub][C];
        int[] remove_causes = new int[N_sub];
        if(check_impossible){
            for(int j = 0; j < C; j++){
                for(int i = 0; i < N; i++){
                    zero_group_matrix[subpop[i]][j] += zero_matrix[i][j];
                }
            }
            for(int j = 0; j < C; j++){
                for(int i = 0; i < N_sub; i++){
                    zero_group_matrix[i][j] = zero_group_matrix[i][j] == 0 ? 0 : 1;
                    remove_causes[i] += 1 - zero_group_matrix[i][j];
                }
            }
        }else{
            for(int i = 0; i < N_sub; i++){
                for(int j = 0; j < C; j++){
                    zero_group_matrix[i][j] = 1;
                }
            }
        }

        // reinitiate after checking impossible
        if(!isAdded){
            // initialize values
            for(int sub  = 0; sub < N_sub; sub++){
                int fix = 0;
                for(int c = 1; c < C; c++){
                    if(zero_group_matrix[sub][c] > 0){
                        fix = c;
                        break;
                    }
                }
                theta_now[sub][fix] = 1;
                double expsum = Math.exp(1.0);
                for(int c = fix + 1; c < C; c++){
                    theta_now[sub][c] = Math.log(rand.nextDouble() * 100.0);
                    expsum += Math.exp(theta_now[sub][c]) * zero_group_matrix[sub][c];
                }
                for(int c = 0; c < C; c++){
                    p_now[sub][c] = Math.exp(theta_now[sub][c]) * zero_group_matrix[sub][c] / expsum;
                }
            }
        }
        // check if specific causes are impossible for all subpopulations
//        int[] zero_pop_vec = new int[C];
//        for(int j = 0; j < C; j++){
//            for(int i = 0; i < N_sub; i++){
//                zero_pop_vec[j] += zero_group_matrix[i][j];
//            }
//            zero_pop_vec[j] = (zero_pop_vec[j] == 0) ? 0 : 1;
//        }


        // first time pnb (naive bayes probability) calculation, note it is inverse of R version
        //        double[][] pnb = new double[N][C];
        double[][] pnb;
        if(!withPhy){
            pnb = insilico.pnb((contains_missing == 1), indic, p_now, subpop, zero_matrix);
        }else{
            int[] g_new = insilico.sampleY(assignment, rand);
            pnb = insilico.pnb((contains_missing == 1), indic, p_now, subpop, g_new, zero_matrix);
        }



        // start loop
        long start = System.currentTimeMillis();
        // popup window under non-unix system
        ProgressPopup popup = new ProgressPopup(this_is_Unix, N_gibbs);

        for(int k = 0; k < N_gibbs; k++){
            if(!this_is_Unix){
                popup.update(k);
            }
            // sample new y vector
            int[] y_new = insilico.sampleY(pnb, rand);
            int[] g_new = insilico.sampleY(assignment, rand);

            // count the appearance of each cause
            int[][] Y = new int[N_sub][C];
            for(int n = 0; n < N; n++){Y[subpop[n]][y_new[n]] += 1;}

            for(int sub = 0; sub < N_sub; sub++){
                // sample mu
                double mu_mean = 0;
                for(int c = 0; c<C; c++) mu_mean += theta_now[sub][c];
                mu_mean = mu_mean / (C-remove_causes[sub]+0.0);
                mu_mean = rngN.nextDouble(mu_mean, Math.sqrt(sigma2_now[sub]/ (C-remove_causes[sub]+0.0)));
                for(int c = 0; c < C; c++){
                    mu_now[sub][c] = mu_mean;
                }

                // sample sigma2
                double shape = (C-remove_causes[sub]-1.0)/2;
                double rate2 = 0;
                for(int c = 0 ; c < C; c++) rate2 += Math.pow(theta_now[sub][c] - mu_now[sub][c] *
                        zero_group_matrix[sub][c], 2);
                sigma2_now[sub] = 1 / rngG.nextDouble( shape, rate2/2 );
                // sample theta
                double[] theta_prev = theta_now[sub];
                theta_now[sub] = insilico.thetaBlockUpdate(jumprange, mu_now[sub],
                        sigma2_now[sub], theta_prev, Y[sub], false, rngN, rand, zero_group_matrix[sub]);
                int j = 0;
                while(j < theta_prev.length){
                    if(theta_now[sub][j] > 0){
                        if(theta_now[sub][j] != theta_prev[j]){
                            naccept[sub] += 1;
                            j = theta_prev.length + 1;
                        }
                    }else{
                        j++;
                    }
                }

                // calculate phat
                double expsum = 0.0;
                for(int c = 0; c < C; c++){
                    expsum += Math.exp(theta_now[sub][c]) * zero_group_matrix[sub][c];
                }
                for(int c = 0; c < C; c++){
                    p_now[sub][c] = Math.exp(theta_now[sub][c]) * zero_group_matrix[sub][c] / expsum;
                }
            }

            if(useProbbase){
                // skip the update step for prob base
            }else{
                insilico.CountCurrent(indic,  y_new);
                if(pool == 0){
                    insilico.TruncBeta_pool(rand, prior_a, prior_b, trunc_min, trunc_max);
                }else if(pool == 1){
                    insilico.TruncBeta(rand, prior_a, prior_b, trunc_min, trunc_max);
                }else if(pool == 2){
                    insilico.TruncBeta2(rand, prior_a, prior_b, trunc_min, trunc_max);
                }
            }

            if(!withPhy){
                pnb = insilico.pnb((contains_missing == 1), indic, p_now, subpop, zero_matrix);
            }else{
                pnb = insilico.pnb((contains_missing == 1), indic, p_now, subpop, g_new, zero_matrix);
            }

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
                System.out.printf("%.2fmin elapsed, %.2fmin remaining \n",
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
                        pnb_mean[d1][d2] += pnb[d1][d2];
                    }
                }
                for(int d1 = 0; d1 < N_sub; d1++){
                    for(int d2 = 0; d2 < C; d2++){
                        p_gibbs[save][d1][d2] = p_now[d1][d2];
                    }
                }
                if(pool == 0){
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
        System.out.println("\nOverall acceptance ratio");
        for(int sub = 0; sub < N_sub; sub++){
            double ratio = naccept[sub] / (N_gibbs+0.0);
            System.out.printf("Sub-population %d : %.4f \n", sub, ratio);
        }


        // The outcome is vector of
        //  0. scalar of N_thin
        //	1. CSMF at each iteration N_sub * C * N_thin
        //  2. Individual prob mean N * C
        //  4. last time configuration N_sub * (C + C + 1)
        //  3. probbase at each iteration

        int N_out = 1 + N_sub * C * N_thin + N * C + N_sub * (C * 2 + 1);
        if(pool == 0){
            N_out +=  N_level * N_thin;
        }else{
            N_out +=  S * C * N_thin;
        }


        // initialize the matrix for output
        double[] parameters = new double[N_out];
        // save CSMF at each iteration
        // counter of which column is
        int counter = 0;
        // save N_thin
        parameters[0] = N_thin;
        counter = 1;
        // save p_gibbs
        for(int sub = 0; sub < N_sub; sub++){
            for(int k = 0; k < N_thin; k++){
                for(int c = 0; c < C; c++){
                    parameters[counter] = p_gibbs[k][sub][c];
                    counter += 1;
                }
            }
        }
        // save pnb_gibbs, need normalize here.
        for(int n = 0; n < N; n++){
            for(int c = 0; c < C; c++){
                parameters[counter] = pnb_mean[n][c]/(N_thin + 0.0);
                counter += 1;
            }
        }
        // save probbase at each iteration
        if(pool != 0){
            for(int c = 0; c < C; c++){
                for(int s = 0; s < S; s++){
                    for(int k = 0; k < N_thin; k++){
                        parameters[counter] = probbase_gibbs[k][s][c];
                        counter += 1;
                    }
                }
            }
        }else{
            for(int k = 0; k < N_thin; k++){
                for(int l = 1; l <= N_level; l++){
                    parameters[counter] = levels_gibbs[k][l-1]; counter+=1;
                }
            }
        }


        // save output without N_thin rows
        // only using the first row
        // save mu_now
        for(int sub = 0; sub < N_sub; sub++){
            for(int c = 0 ; c < C; c++){
                parameters[counter] = mu_now[sub][c];
                counter += 1;
            }
        }
        // save sigma2_now
        for(int sub = 0; sub < N_sub; sub++){
            parameters[counter] = sigma2_now[sub];
            counter += 1;
        }
        // save theta_now
        for(int sub = 0; sub < N_sub; sub++){
            for(int c = 0 ; c < C; c++){
                parameters[counter] = theta_now[sub][c];
                counter += 1;
            }
        }

//	else{
//			// save output without N_thin rows
//			// only using the last row
//			// starting from the first column
//			counter = 0;
//			for(int sub = 0; sub < N_sub; sub++){
//				for(int c = 0 ; c < C; c++){
//					parameters[N_out_1 - 1][counter] = mu_now[sub][c];
//					counter += 1;
//				}
//			}
//			// save sigma2_now
//			for(int sub = 0; sub < N_sub; sub++){
//				parameters[N_out_1 - 1][counter] = sigma2_now[sub];
//				counter += 1;
//			}
//			// save theta_now
//			for(int sub = 0; sub < N_sub; sub++){
//				for(int c = 0 ; c < C; c++){
//					parameters[N_out_1 - 1][counter] = theta_now[sub][c];
//					counter += 1;
//				}
//			}
//		}
        System.out.println("Organizing output, might take a moment...");

//		return;
        return(parameters);
    }

    /**
     *
     * @param DontAsk: M by 8 matrix (java index + 1, 0 if not exist)
     * @param AskIf: M by 2 matrix (java index + 1, 0 if not exist)
     * @param data: N by M matrix
     * @param warningfile: warning file
     * @return
     */
    public static double[][] Datacheck(int[][] DontAsk, int[][] AskIf, int[] zero_to_missing_list,  double[][] data, String[] ids, String[] symps, String warningfile)
            throws FileNotFoundException, UnsupportedEncodingException {

        boolean is_this_missing;
        boolean write = ids.length > 1;
        boolean[] set_to_missing = new boolean[data[0].length];
        if(zero_to_missing_list[0] != 0){
            for(int i : zero_to_missing_list){
                set_to_missing[i-1] = true;
            }
        }

        PrintWriter writer = new PrintWriter(warningfile, "UTF-8");
        if(write){
            writer.println("Warning log built for InSilicoVA");
        }
        for(int i = 0; i < data.length; i++){
           for(int k = 0; k < 2; k++){
               for(int j = 0; j < data[i].length; j++ ){
                   if(data[i][j] != -1){
                        for(int t = 0; t < DontAsk[j].length; t++){
                            if(DontAsk[j][t] == 0){continue;}
                            if(data[i][DontAsk[j][t] - 1] == 1){
                                /** set to missing only if (a) Dont ask is an demographic indicator, and (b) this symptom is not an demographic indicator **/
                                is_this_missing = set_to_missing[DontAsk[j][t] - 1] ;//& (! set_to_missing[j]);
                                data[i][j] = is_this_missing ? -1 : 0;
                                if(write) writer.println(ids[i] + " " + symps[j]+ " value inconsistent with " +  symps[DontAsk[j][t]] + " - cleared or set to missing in working file");
                                // break if specified to set to missing
                                if(is_this_missing) break;
                            }
                        }
                       for(int t = 0; t < AskIf[j].length; t++){
                           if(AskIf[j][t] == 0 | data[i][j] != 1){continue;}
                           data[i][AskIf[j][t] - 1] = 1;
                           if(write) writer.println(ids[i] + " " + symps[j]+ " not flagged in category " + symps[DontAsk[j][t]] + " - updated in working file");
                       }
                   }else if(data[i][j] == -1){
                       for(int t = 0; t < DontAsk[j].length; t++){
                           if(DontAsk[j][t] == 0){continue;}
                           if(data[i][DontAsk[j][t] - 1] == 1){
                               is_this_missing = set_to_missing[DontAsk[j][t] - 1] ;//& (! set_to_missing[j]);
                               data[i][j] = is_this_missing ? -1 : 0;
                               if(write) writer.println(ids[i] + " " + symps[j]+ " value inconsistent with " + symps[DontAsk[j][t]] + " - set to No or Missing in working file");
                               // break if specified to set to missing
                               if(is_this_missing) break;
                           }
                       }
                   }
               }
           }
        }
        writer.close();
        return(data);
    }

    /**
     *
     * @param DontAsk: M by 8 matrix (java index + 1, 0 if not exist)
     * @param AskIf: M by 2 matrix (java index + 1, 0 if not exist)
     * @param data: N by M matrix
     * @return
     */
    public static double[][] Datacheck(int[][] DontAsk, int[][] AskIf, int[] zero_to_missing_list, double[][] data) throws FileNotFoundException, UnsupportedEncodingException {
        String[] ids = new String[0];
        String[] symps = new String[0];
        String warningfile = "warnings.txt";
        return(Datacheck(DontAsk, AskIf, zero_to_missing_list, data, ids, symps, warningfile));
    }

    public static double[][] IndivProb(double[][] data, int[][] impossible, double[]csmf0, int[] subpop,
                                       double[] condprob0, double p0, double p1,
                                       int Nsub, int Nitr, int C, int S) {
        int N = data.length;

        /**   --------------- to bypass bug in rJava package  ----------------------- **/
        // reconstruct csmf and condprob csmf[i, j, k] <-> csmf0[k*d1*d2 + j*d1 + i]
        double[][][] csmf = new double[Nsub][Nitr][C];
        double[][][] condprob = new double[Nitr][S][C];
        for(int i = 0; i < Nsub; i++){
            for(int j = 0; j < Nitr; j++){
                for(int k = 0; k < C; k++){
                    csmf[i][j][k] = csmf0[k * Nsub * Nitr + j * Nsub + i];
                }
            }
        }
        for(int i = 0; i < Nitr; i++){
            for(int j = 0; j < S; j++){
                for(int k = 0; k < C; k++){
                    condprob[i][j][k] = condprob0[k * Nitr * S + j * Nitr + i];
                }
            }
        }
        /**   --------------- reconstruction done ----------------------- **/

        int[][] zero_matrix = new int[N][C];
        for(int i = 0; i < N; i++){
            for(int j = 0; j < C; j++){
                zero_matrix[i][j] = 1;
            }
        }
        if(impossible[0].length == 2){
                for(int i = 0; i < N; i++){
                    for(int k = 0; k < impossible.length; k++){
                        if(data[i][impossible[k][1] - 1] == 1){
                            zero_matrix[i][impossible[k][0] - 1] = 0;
                        }
                    }
                }
            }
        double[][][] allp = new double [N][C][Nitr];

        for(int i = 0; i < N; i++){
            // get subpop index for death i
            int sub = subpop[i] - 1;
            for(int t = 0; t < Nitr; t++) {

                /** ------------------------------------------------**/
                double sum = 0;
                for (int c = 0; c < C; c++) {
                    allp[i][c][t] = csmf[sub][t][c] * zero_matrix[i][c];
                    for (int s = 0; s < S; s++) {
                        if(condprob[t][s][c] == 0){System.out.print(".");}

                        if (data[i][s] == 1) {
                            allp[i][c][t] *= condprob[t][s][c];
                        } else if (data[i][s] == 0) {
                            allp[i][c][t] *= 1 - condprob[t][s][c];
                        }
                    }
                    sum += allp[i][c][t];
                }
                for (int c = 0; c < C; c++) {
                    allp[i][c][t] = (sum == 0) ? allp[i][c][t] : allp[i][c][t] / sum;
                }
                /** ------------------------------------------------**/
            }
        }


        double[][] out = new double[N*4][C];
        for(int i = 0; i < N; i++){
            for(int c = 0; c < C; c++){
                out[i][c] = MathUtil.getMean(allp[i][c]);
                out[i + N][c] = MathUtil.getPercentile(allp[i][c], 0.5);
                out[i + N * 2][c] = MathUtil.getPercentile(allp[i][c], p0);
                out[i + N * 3][c] = MathUtil.getPercentile(allp[i][c], p1);
            }
        }
        return(out);
    }

    public static double[][] AggIndivProb(double[][] data, int[][] impossible, double[]csmf0, int[] subpop,
                                          double[] condprob0, int[] group, int Ngroup, double p0, double p1,
                                          int Nsub, int Nitr, int C, int S) {
        int N = data.length;

        /**   --------------- to bypass bug in rJava package  ----------------------- **/
        // reconstruct csmf and condprob csmf[i, j, k] <-> csmf0[k*d1*d2 + j*d1 + i]
        double[][][] csmf = new double[Nsub][Nitr][C];
        double[][][] condprob = new double[Nitr][S][C];
        for(int i = 0; i < Nsub; i++){
            for(int j = 0; j < Nitr; j++){
                for(int k = 0; k < C; k++){
                    csmf[i][j][k] = csmf0[k * Nsub * Nitr + j * Nsub + i];
                }
            }
        }
        for(int i = 0; i < Nitr; i++){
            for(int j = 0; j < S; j++){
                for(int k = 0; k < C; k++){
                    condprob[i][j][k] = condprob0[k * Nitr * S + j * Nitr + i];
                }
            }
        }
        /**   --------------- reconstruction done ----------------------- **/

        int[][] zero_matrix = new int[N][C];
        for(int i = 0; i < N; i++){
            for(int j = 0; j < C; j++){
                zero_matrix[i][j] = 1;
            }
        }
        if(impossible[0].length == 2){
            for(int i = 0; i < N; i++){
                for(int k = 0; k < impossible.length; k++){
                    if(data[i][impossible[k][1] - 1] == 1){
                        zero_matrix[i][impossible[k][0] - 1] = 0;
                    }
                }
            }
        }
        double[][][] allp = new double [Ngroup][C][Nitr];
        int[] size = new int[Ngroup];

        for(int i = 0; i < N; i++){
            // get subpop index for death i
            int sub = subpop[i] - 1;
            int g = group[i] - 1;
            size[g] ++;

            for(int t = 0; t < Nitr; t++) {
                double[] tmp = new double[C];
                /** ------------------------------------------------**/
                double sum = 0;
                for (int c = 0; c < C; c++) {
                    tmp[c] = csmf[sub][t][c] * zero_matrix[i][c];
                    for (int s = 0; s < S; s++) {
                        if (data[i][s] == 1) {
                            tmp[c] *= condprob[t][s][c];
                        } else if (data[i][s] == 0) {
                            tmp[c] *= 1 - condprob[t][s][c];
                        }
                    }
                    sum += tmp[c];
                }
                for (int c = 0; c < C; c++) {
                    allp[g][c][t] += (sum == 0) ? tmp[c] : tmp[c] / sum;
                }
                /** ------------------------------------------------**/
            }
        }


        double[][] out = new double[Ngroup*4][C+1];
        for(int i = 0; i < Ngroup; i++){
            for(int c = 0; c < C; c++){
                out[i][c] = MathUtil.getMean(allp[i][c]) / (size[i] + 0.0);
                out[i + Ngroup][c] = MathUtil.getPercentile(allp[i][c], 0.5) / (size[i] + 0.0);
                out[i + Ngroup * 2][c] = MathUtil.getPercentile(allp[i][c], p0) / (size[i] + 0.0);
                out[i + Ngroup * 3][c] = MathUtil.getPercentile(allp[i][c], p1) / (size[i] + 0.0);
            }
            out[i][C] = size[i];
            out[i + Ngroup][C] = size[i];
            out[i + Ngroup * 2][C] = size[i];
            out[i + Ngroup * 3][C] = size[i];
        }
        return(out);
    }

}