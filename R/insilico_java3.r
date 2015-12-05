## Created on Oct 1, 2015


#' Implement InSilicoVA methods with more flexible customization
#' 
#' This function implements InSilicoVA model. This is the developer's version of 
#' InSilicoVA with more flexibility in customized input. 
#' The InSilicoVA model is fitted
#' with MCMC implemented in Java. For more detail, see the paper on
#' \url{http://arxiv.org/abs/1411.3042}.
#' 
#' For Windows user, this function will produce a popup window showing the
#' progress. For Mac and Unix user, this function will print progress messages
#' on the console. Special notice for users using default R GUI for mac, the
#' output will not be printed on console while the function is running, and
#' will only be printed out after it is completed. Thus if you use a Mac, we
#' suggest using either RStudio for mac, or running R from terminal.
#' 
#' The chains could be set to run automatically longer. If set
#' \code{auto.length} to be TRUE, the chain will assess convergence after
#' finishing the length K chain input by user using Heidelberger and Welch's
#' convergence diagnostic. If convergence is not reached, the chain will run
#' another K iterations and use the first K iterations as burn-in. If the chain
#' is still not converged after 2K iterations, it will proceed to another 2K
#' iterations and again use the first 2K iterations as burn-in. If convergence
#' is still not reached by the end, it will not double the length again to
#' avoid heavy memory use. A warning will be given in that case. The extended
#' chains will be thinned in the same way.
#' 
#' For more detail of model specification, see the paper on
#' \url{http://arxiv.org/abs/1411.3042}.
#' 
#' 
#' @param data The original data to be used. It is suggested to use similar
#' input as InterVA4, with the first column being death IDs. The only
#' difference in input is InsilicoVA takes three levels: ``present'',
#' ``absent'', and ``missing (no data)''. Similar to InterVA software,
#' ``present'' symptoms takes value ``Y''; ``absent'' symptoms take take value
#' ``NA'' or ``''. For missing symptoms, e.g., questions not asked or answered
#' in the original interview, corrupted data, etc., the input should be coded
#' by ``.'' to distinguish from ``absent'' category. The order of the columns does
#' not matter as long as the column names are correct. It can also include more 
#' unused columns than the standard InterVA4 input. But the first column should be 
#' the death ID.
#' @param isNumeric Indicator if the input is already in numeric form. If the
#' input is coded numerically such that 1 for ``present'', 0 for ``absent'',
#' and -1 for ``missing'', this indicator could be set to True to avoid
#' conversion to standard InterVA format.
#' @param useProbbase Logical indicator. If TRUE, then use InterVA conditional
#' probability table without re-estimating.
#' @param keepProbbase.level Logical indicator when \code{useProbbase} is
#' FALSE. If TRUE, then only estimate the InterVA's conditional probability
#' table interpretation table; if FALSE, estimate the whole conditional
#' probability matrix. Default to TRUE.
#' @param cond.prob.touse Customized conditional probability table to use.
#' Currently only accepting the same configuration as InterVA-4 software. It
#' should be a matrix of 245 rows of symptoms and 60 columns of causes,
#' arranged in the same order as in InterVA-4 specification. The elements in
#' the matrix should be the conditional probability of corresponding symptom
#' given the corresponding cause, represented in alphabetic form indicating
#' levels. For example input, see \code{\link{condprob}}
#' @param datacheck Logical indicator for whether to check the data satisfying
#' InterVA rules. Default set to be TRUE. If \code{warning.write} is set to
#' true, the inconsistent input will be logged in file warnings.txt. It's
#' strongly suggested to be set to TRUE.
#' @param warning.write Logical indicator for whether to save the changes made
#' to data input by \code{datacheck}. If set to TRUE, the changes will be
#' logged in file warnings.txt in current working directory.
#' @param external.sep Logical indicator for whether to separate out external
#' causes first. Default set to be TRUE. If set to TRUE, the algorithm will
#' estimate external causes, e.g., traffic accident, accidental fall, suicide,
#' etc., by checking the corresponding indicator only without considering other
#' medical symptoms. It is strongly suggested to set to be TRUE.
#' @param length.sim Number of iterations to run. Default to be 4000.
#' @param thin Proportion of thinning for storing parameters. For example, if
#' thin = k, the output parameters will only be saved every k iterations.
#' Default to be 10
#' @param burnin Number of iterations as burn-in period. Parameters sampled in
#' burn-in period will not be saved.
#' @param auto.length Logical indicator of whether to automatically increase
#' chain length if convergence not reached.
#' @param conv.csmf Minimum CSMF value to check for convergence if auto.length
#' is set to TRUE. For example, under the default value 0.02, all causes with
#' mean CSMF at least 0.02 will be checked for convergence.
#' @param jump.scale The scale of Metropolis proposal in the Normal model.
#' Default to be 0.1.
#' @param levels.prior Vector of prior expectation of conditional probability
#' levels. They do not have to be scaled. The algorithm internally calibrate
#' the scale to the working scale through \code{levels.strength}. If NULL the
#' algorithm will use InterVA table as prior.
#' @param levels.strength Scaling factor for the strength of prior beliefs in
#' the conditional probability levels. Larger value constrain the posterior
#' estimates to be closer to prior expectation. Defult value 1 scales
#' \code{levels.prior} to a suggested scale that works empirically.
#' @param trunc.min Minimum possible value for estimated conditional
#' probability table. Default to be 0.0001
#' @param trunc.max Maximum possible value for estimated conditional
#' probability table. Default to be 0.9999
#' @param subpop This could be the column name of the variable in data that is to
#' be used as sub-population indicator, or a list of column names if more than one 
#' variable are to be used. Or it could be a vector of sub-population assignments 
#' of the same length of death records. It could be numerical indicators or character 
#' vectors of names. 
#' @param java_option Option to initialize java JVM. Default to ``-Xmx1g'',
#' which sets the maximum heap size to be 1GB. If R produces
#' ``java.lang.OutOfMemoryError: Java heap space'' error message, consider
#' increasing heap size using this option, or one of the following: (1)
#' decreasing \code{length.sim}, (2) increasing \code{thin}, or (3) disabling
#' \code{auto.length}.
#' @param seed Seed used for initializing sampler. The algorithm will produce
#' the same outcome with the same seed in each machine.
#' @param phy.code A matrix of physician assigned cause distribution. The
#' physician assigned causes need not be the same as the list of causes used in
#' InSilicoVA and InterVA-4. The cause list used could be a higher level
#' aggregation of the InSilicoVA causes. See \code{phy.cat} for more detail.
#' The first column of \code{phy.code} should be death ID that could be matched
#' to the symptom dataset, the following columns are the probabilities of each
#' cause category used by physicians.
#' @param phy.cat A two column matrix describing the correspondence between
#' InSilicoVA causes and the physician assigned causes. Note each InSilicoVA
#' cause (see \code{causetext}) could only correspond to one physician assigned
#' cause. See \code{SampleCategory} for an example. 'Unknown' category should
#' not be included in this matrix.
#' @param phy.unknown The name of the physician assigned cause that correspond
#' to unknown COD.
#' @param phy.external The name of the physician assigned cause that correspond
#' to external causes. This will only be used if \code{external.sep} is set to
#' TRUE. In that case, all external causes should be grouped together, as they
#' are assigned deterministically by the corresponding symptoms.
#' @param phy.debias Fitted object from physician coding debias function (see
#' \code{physician_debias}) that overwrites all physician coding parameters.
#' @param dev.customization default to be FALSE
#' @param Probbase_by_symp.dev default to be FALSE
#' @param probbase.dev default to be NULL
#' @param table.dev default to be NULL
#' @param gstable.dev default to be NULL
#' @param nlevel.dev default to be NULL
#' @param prob.order.dev default to be NULL
#' @return \item{id}{A vector of death ID. Note the order of the ID is in
#' general different from the input file. See \code{report} for organizing the
#' report.}
#' 
#' \item{indiv.prob}{Matrix of individual mean cause of death distribution.
#' Each row corresponds to one death with the corresponding ID.}
#' 
#' \item{csmf}{Matrix of CSMF vector at each iterations after burn-in and
#' thinning. Each column corresponds to one cause.}
#' 
#' \item{conditional.probs}{If the model is estimated with
#' \code{keepProbbase.level} = TRUE, this value gives a matrix of each
#' conditional probability at each level at each iterations. Each column
#' corresponds to one level of probability. If \code{keepProbbase.level} =
#' FALSE, this value gives a three-dimensional array. If \code{UseProbbase} =
#' TRUE, the value will be set to NULL. See \code{report} for more analysis.}
#' 
#' \item{missing.symptoms}{Vector of symptoms missing from all input data.}
#' 
#' \item{external}{Logical indicator of whether the model is fitted with
#' external causes separated calculated.}
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{plot.insilico}}, \code{\link{summary.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark(2014) \emph{Probabilistic
#' cause-of-death assignment using verbal autopsies},
#' \url{http://arxiv.org/abs/1411.3042} \cr \emph{Working paper no. 147, Center
#' for Statistics and the Social Sciences, University of Washington}
#' @keywords InSilicoVA
#' @examples
#' 
#' \dontrun{
#' # load sample data together with sub-population list
#' data(RandomVA1)
#' # extract InterVA style input data
#' data <- RandomVA1$data
#' # extract sub-population information. 
#' # The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
#' subpop <- RandomVA1$subpop
#' 
#' # run without subpopulation
#' fit1<- insilico( data, subpop = NULL,  
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE)
#' 
#' # re-run with subpopulation
#' fit2<- insilico( data, subpop = subpop, 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE)
#'  
#' # note a different ways to specify subpopulation
#' HIV_indicator <- RandomVA1$subpop
#' data_more <- cbind(data, HIV_indicator)
#' fit2<- insilico( data_more, subpop = "HIV_indicator", 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE)
#'
#' # When more than one variables are used to divide the population
#' # for example, if we shuffle the HIV indicator to create another fake grouping
#' # and we want to use the combination of both as a subpopulation indicator
#' HIV_indicator <- RandomVA1$subpop
#' HIV_indicator2 <- sample(RandomVA1$subpop)
#' data_more <- cbind(data, HIV_indicator, HIV_indicator2)
#' fit2<- insilico( data_more, subpop = list("HIV_indicator", "HIV_indicator2"), 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE)
#'
#'
#' # run without re-sampling conditional probabilities
#' fit3<- insilico( data, subpop = subpop, 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, useProbbase = TRUE)
#' }
#' 
#' @export insilico
insilico.dev <- function(data, isNumeric = FALSE,useProbbase = FALSE, keepProbbase.level = TRUE,  cond.prob.touse = NULL,datacheck = TRUE, warning.write = FALSE, external.sep = TRUE, length.sim = 4000, thin = 10, burnin = 2000, auto.length = TRUE, conv.csmf = 0.02, jump.scale = 0.1, levels.prior = NULL, levels.strength = 1, trunc.min = 0.0001, trunc.max = 0.9999, subpop = NULL, java_option = "-Xmx1g", seed = 1, phy.code = NULL, phy.cat = NULL, phy.unknown = NULL, phy.external = NULL, phy.debias = NULL, dev.customization = FALSE, Probbase_by_symp.dev = FALSE, probbase.dev = NULL, table.dev = NULL, gstable.dev = NULL, nlevel.dev = NULL, prob.order.dev = NULL){ 
	
#############################################################################
#############################################################################
## InSilico VA -  helper functions 
##
## author: Richard Li 
## date: 05/11/2014
#############################################################################
InterVA.table <- function(min){
###########################################################
# function to return the interVA conversion table for alphabetic levels
# also change the smallest value from 0 to user input		
# @param:
#		min: minimum level
# @values:
# 		vector of interVA4 levels
	return(c(1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 
			  0.001, 0.0005, 0.0001, 0.00001, min))
}	
scale.vec <- function(aaa, scale = NULL, scale.max = NULL, reverse = TRUE){
###########################################################
# function to change the scale of a given vector
# Note: need to reverse the order! so A has larger value, F has smaller value
# @param:
#	    aaa       : vector to be scaled
# 		scale     : sum of the vector after scaling
#       scale.max : max of the vector after scaling
#		reverse   : whether the order should be reversed
# @values:	
#  		scaled vector
	if(reverse){aaa <- max(aaa) + 1 - aaa}
	bbb <- aaa
	if(!is.null(scale)) return(bbb/sum(bbb) * scale)
	if(!is.null(scale.max)) return(bbb * scale.max / max(bbb))
}

scale.vec.inter <- function(aaa, scale = NULL, scale.max = NULL){
###########################################################
# function to map ordered input vector to InterVA alphabetic scale
# Note:  to avoid 0, change the last one to 0.000001 here
# @param:
#	    aaa       : vector to be scaled
# 		scale     : sum of the vector after scaling
#       scale.max : max of the vector after scaling
#		reverse   : whether the order should be reversed
# @values:	
#  
	dist <- InterVA.table(0.000001)
	if(length(aaa) != length(dist)){stop("dimension of probbase prior not correct")}
	bbb <- dist[order(aaa)]
	if(!is.null(scale)) return(bbb/sum(bbb) * scale)
	if(!is.null(scale.max)) return(bbb * scale.max / max(bbb))
}

change.inter <- function(x, order = FALSE){
# function to translate alphebatic matrix into numeric matrix or order matrix
# @param:
# 	x      : alphabetic matrix 
#	order  : whether to change the matrix into order matrix
# @values:
#	numeric matrix by InterVA probbase rules, or the order matrix
	a <- dim(x)[1]
	b <- dim(x)[2]
	if(is.null(a)){
		y <- rep(0,length(x))
	}else{
		y <- matrix(0, a, b)
	}  	
	inter.table <- InterVA.table(0)
	y[x == "I"] <- 1
    y[x == "A+"] <- 0.8
    y[x == "A"] <- 0.5
    y[x == "A-"] <- 0.2
    y[x == "B+"] <- 0.1
    y[x == "B"] <- 0.05
    y[x == "B-"] <- 0.02
    y[x == "B -"] <- 0.02
    y[x == "C+"] <- 0.01
    y[x == "C"] <- 0.005
    y[x == "C-"] <- 0.002
    y[x == "D+"] <- 0.001
    y[x == "D"] <- 5e-4
    y[x == "D-"] <- 1e-4
    y[x == "E"] <- 1e-5
    y[x == "N"] <- 0
    y[x == ""] <- 0

    if(order){
    	for(i in 1: length(inter.table)){
			y[y == inter.table[i] ] <- i
		}
    }
    if(!is.null(a)){
		y <- matrix(y, a, b)
	}  	
    return(y)   
}

cond.initiate <- function(probbase.order, expIni, Inter.ini, min, max){
# Randomly initialize probbase from order matrix
# @param:
# 	probbase.order : order matrix
#   expIni         : initialize to exponential of uniform
#   Inter.ini      : initialize to InterVA probbase values
#   min            : minimum of probbase values
#	max            : maximum of probbase values
# @values:
#	new probbase matrix

		# move initial values away from 0 and 1
		if(min == 0) min <- 0.01
		if(max == 1) max <- 0.99		
		# find the number of levels
		nlevel <- max(probbase.order)
		# if the random levels are in log-linear fashion
		if(expIni == TRUE){
		randomlevels <- sort(exp(runif(nlevel, log(min), log(max))), decreasing = TRUE)	
		}else{
		randomlevels <- sort(runif(nlevel, min, max), decreasing = TRUE)		  
		}
		# if the random levels are initialized proportional to InterVA intervals
		if(Inter.ini){
			randomlevels <- InterVA.table(0)
			randomlevels <- (randomlevels * (max - min)) + min
		}
		# change order matrix into actual values
		probrandom = probbase.order
		for(i in 1:nlevel){
			probrandom[which(probbase.order == i)] <- randomlevels[i]
		}

		return(probrandom)	
}

		
removeBad <- function(data, is.numeric, subpop){
###########################################################
# function to remove data with no sex/age indicators 
# @param:
#		data        : as in main function
#		is.numeric  : as in main function
#		subpop      : as in main function
# @values:
# 		a list of data and subpop after removing bad data
	if(!is.numeric){
		data.num <- matrix(0, dim(data)[1], dim(data)[2])		
		for(i in 2:dim(data)[2]){
			temp <- data[,i]
			temp.ind <- which(temp == "Y")
			data.num[temp.ind, i] <- 1
		}
	}else{
		data.num <- data
	}
	err <- NULL
	for(i in 1:dim(data.num)[1]){
		if(sum(data.num[i, 2:8]) < 1 ||
			sum(data.num[i, 9:10]) < 1 ||
			sum(data.num[i, 23:223]) < 1){
			err <- c(err, i)
		}
	}
	if(is.null(err)) return(list(data = data, subpop = subpop))
	if(!is.null(subpop)) return(list(data = data[-err, ], 
									 subpop = subpop[-err]))
	if(is.null(subpop)) return(list(data = data[-err, ], subpop = NULL))
}

datacheck.interVA <- function(id, indic, missing.all, external.sep, warning.write){
###########################################################
# function to perform data check as in InterVA4 
# @param:
#		id             : vector of IDs
#		indic          : matrix of indicators
#		missing.all    : as in main function
#		external.sep   : as in main function
#		warning.write  : as in main function
# @values:
# 		revised indic matrix	
		if(warning.write){
			cat(paste("Warning log built for InterVA", Sys.time(), "\n"),file="warnings.txt",append = FALSE) 
		}
		data("probbase", envir = environment())
		probbase <- get("probbase", envir  = environment())
		data("causetext", envir = environment())
		causetext <- get("causetext", envir  = environment())

	    probbase <- as.matrix(probbase)
	    Input <- as.matrix(indic)
	    symps <- probbase[2:246, 2]
	    if(external.sep){
	    	symps <- symps[-c(211:222)]
	    }
	    if(!is.null(missing.all)){
	    	symps <- symps[-missing.all]
	    }
	    colnames(Input) <- symps
	    for(i in 1:dim(Input)[1]){
	        input.current <- as.numeric(Input[i,])
	        ## Repeat twice the check of "ask if" and "don't ask".
	        ## If there is contradictory with "ask if" or "don't ask", follow the following rules:
	        ## If B is the "don't ask" for A but B has value 1 --> make sure A has value 0;
	        ## If B is the "ask if" for A but B has value 0 --> change B into value 1
	        for(k in 1:2){
	            for(j in 1:(dim(Input)[2] )){
	                if(input.current[j] == 1 ){
	                    Dont.ask <- probbase[j + 1, 4:11]
	                    Ask.if <- probbase[j + 1, 12]
	                    Dont.ask.list <- input.current[match(toupper(Dont.ask), toupper(colnames(Input)))]
	                    Dont.ask.list[ is.na(Dont.ask.list)] <- 0
	                    
	                    if( max( Dont.ask.list ) > 0 ){
	                    	input.current[j] <- 0
	                    	if(warning.write){
	                    		cat(id[i], "   ", paste(probbase[j+1, 2], "  value inconsistent with ", Dont.ask[which(Dont.ask.list > 0)], " - cleared in working file \n"), file="warnings.txt", append=TRUE)
	                    	}
	                    }
	                    if( !is.na(match(Ask.if, colnames(Input)))  ){
	                        if(input.current[match(Ask.if, colnames(Input) )] < 1){
	                            input.current[match(Ask.if, colnames(Input) )] <- 1
	                             if(warning.write){
                            	    cat(id[i], "   ", paste(probbase[j+1, 2], "  not flagged in category ", Ask.if, " - updated in working file \n"), file="warnings.txt", append=TRUE)
                           		 }
	                        }
	                    }
	                }
	            }
	        }
	   Input[i, ] <- input.current          
	}
	Input
}

impute <- function(indic.w.missing, y.new, cond.prob, S, subbelong, missing.imp){
###########################################################
# function to impute missing data (NOT USED)
# @param:
#		indic.w.missing : matrix of indicators
#		y.new           : current cause vector
#		cond.prob       : cond prob matrix
#		S 				: number of symptoms
#		subbelong       : vector of sub-population membership
#		missing.imp     : as in the main function	
# @values:
# 		imputed indic matrix		
	## if no imputing
	if(missing.imp == "none"){
		return(indic.w.missing)

	## if impute for all missings
	}else if(missing.imp == "all"){
		## arrange indic into vectors, one death after another
		indic.w.missing <- t(indic.w.missing)
		## find where are the missing values
		where.miss <- which(indic.w.missing < 0)
	
	# if impute for all missings except sub-population specific total missing
	}else if(missing.imp == "sub"){
		## arrange indic into vectors, one death after another
		indic.w.missing <- t(indic.w.missing)
		## find where are the missing values
		where.miss <- which(indic.w.missing == -1)
	}

	if(length(where.miss)  == 0) return(indic.w.missing)
	## find the causes for each of the missing values
	what.cause <- y.new[trunc((where.miss-0.1)/S) + 1]
	## find the corresponding cond.prob
	getprob <- function(i, where.miss, what.cause, cond.prob, S){
		which.symp <- where.miss[i] %% S
		if(which.symp == 0) which.symp <- S
		cond.prob[which.symp, what.cause[i]]
	}
	prob <- sapply(seq(1:length(where.miss)), getprob, 
					where.miss, what.cause, cond.prob, S)
	## sample new indic.vec
	samp <- rbinom(n = length(where.miss), size = 1, 
					prob = prob)
	indic.w.missing[where.miss] <- samp
	return(t(indic.w.missing))
}
  
removeExt <- function(data, prob.orig, is.Numeric, subpop, subpop_order_list, external.causes, external.symps){
###########################################################
# function to remove external causes/symps and assign deterministic deaths
# @param:
#	   data
#      prob.orig: directly from InterVA4
# @values:
#	   data: after removing external death
#	   prob.orig: after deleting external symptoms
#	   exts: external death list		
	extSymps <- external.symps
	extCauses <- external.causes
	# extract subset of data for external symptoms
	N.all <- dim(data)[1]
	extData <- data[, extSymps + 1]
	if(is.Numeric){
		neg <- 0
		pos <- 1
	}else{
		neg <- ""
		pos <- "Y"
	}
	ext.where <- which(apply(extData, 1, function(x){
									length(which(x == pos)) }) > 0)
	extData <- as.matrix(extData[ext.where, ])
	ext.id <- data[ext.where, 1]
	ext.sub <- subpop[ext.where]
	# a smaller scale datacheck for external causes
	extData[which(extData[,3] == pos), 2] <- neg
	extData[which(extData[,4] == pos), c(2,3,5,6)] <- neg
	extData[which(extData[,5] == pos), c(2,3,4)] <- neg
	extData[which(extData[,6] == pos), c(8,10)] <- neg
	extData[which(extData[,8] == pos), c(2,3,4,5)] <- neg
	extData[which(extData[,9] == pos), c(2,3,4,8,10)] <- neg
	extData[which(extData[,10] == pos),c(2,4,5,6,8)] <- neg
	extData[which(extData[,11] == pos),c(3,4,8,9)] <- neg
	extData[which(extData[,12] == pos),c(4,7,8)] <- neg
	
	# initialize with all "unspecified ext causes"
	ext.cod <- rep(extCauses[11], length(ext.id))
	# begin checking symptoms
	# road traffic
	ext.cod[which(extData[,2] == pos)] <- extCauses[1]
	# non-road transport
	ext.cod[which(extData[,3] == pos)] <- extCauses[2]
	# accident fall
	ext.cod[which(extData[,4] == pos)] <- extCauses[3]
	# drowning
	ext.cod[which(extData[,5] == pos)] <- extCauses[4]
	# burns
	ext.cod[which(extData[,6] == pos)] <- extCauses[5]
	# assault
	ext.cod[which(extData[,7] == pos)] <- extCauses[10]
	# animal bite
	ext.cod[which(extData[,8] == pos)] <- extCauses[6]
	# force of nature
	ext.cod[which(extData[,9] == pos)] <- extCauses[7]
	# poison
	ext.cod[which(extData[,10] == pos)] <- extCauses[8]
	# inflict
	ext.cod[which(extData[,11] == pos)] <- extCauses[10]
	# suicide
	ext.cod[which(extData[,12] == pos)] <- extCauses[9]

	# delete death confirmed external
	data <- data[-ext.where, ]

	data <- data[, -(extSymps + 1)]
	# delete the causes from probbase
	prob.orig <- prob.orig[ -(extSymps), -(extCauses)]
	if(!is.null(subpop)){
		ext.csmf <- vector("list", length(subpop_order_list))
		for(i in 1:length(ext.csmf)){
			ext.csmf[[i]] <- rep(0, length(extCauses))
			ext.cod.temp <- ext.cod[which(ext.sub == subpop_order_list[i])]
			if(!is.null(ext.cod.temp)){
				for(j in 1:length(extCauses)){
					ext.csmf[[i]][j] <- length(which(ext.cod.temp == extCauses[j]))
				}
				ext.csmf[[i]] <- ext.csmf[[i]]/length(which(subpop == subpop_order_list[i]))		
			}
		}
	}else{
		ext.csmf <- rep(0, length(extCauses))
		for(i in 1:length(extCauses)){
			ext.csmf[i] <- length(which(ext.cod == extCauses[i]))
		}
		ext.csmf <- ext.csmf/N.all		
	}
	return(list(data = data, 
				subpop = subpop[-ext.where],
				prob.orig = prob.orig, 
				ext.sub  = ext.sub,
				ext.id = ext.id, 
				ext.cod = ext.cod,
				ext.csmf = ext.csmf))
}

ParseResult <- function(N_sub.j, C.j, S.j, N_level.j, pool.j, fit){
###########################################################
# function to parse results from Java into correct place
# @param:
#	   various java arguments
# 	   Java output
# @values:
#		list of variables parsed   
   	counter <- 1
    csmf.sub <- NULL
    p.hat <- NULL
    probbase.gibbs <- NULL
    levels.gibbs <- NULL

    # extract N_thin
    N_thin = fit[1]
    counter = counter + 1

    # extract CSMF
    if(N_sub.j > 1){
        csmf.sub <- vector("list", N_sub.j)
        for(sub in 1:N_sub.j){
            csmf.sub[[sub]] <- fit[counter : (counter + C.j * N_thin - 1)]
            csmf.sub[[sub]] <- matrix(csmf.sub[[sub]], nrow = N_thin, ncol = C.j, 
            						  byrow = TRUE)
            counter <- counter + C.j * N_thin
        }
    }else{
        p.hat <- fit[counter : (counter + C.j * N_thin- 1)]
        p.hat <- matrix(p.hat, nrow = N_thin, ncol = C.j, byrow = TRUE)
        counter <- counter + C.j * N_thin
    }
    
    # extract individual probabilities
    p.indiv <- fit[counter : (counter + N.j * C.j - 1)]
    p.indiv <- matrix(p.indiv, nrow = N.j, ncol = C.j, byrow = TRUE)
    counter <- counter + N.j * C.j

    # extract probbase 
    if(pool.j != 0){
        probbase.gibbs <- fit[counter:(counter + S.j * C.j * N_thin - 1)]
        # array(..., dim =c(a,b,c))
        # what it does it for each c, fill by column
        # i.e. c(x[1,1,1], x[2, 1, 1], x[1, 2, 1], ...) 
        # i.e. in Java, loop in the order of C.j -> S.j -> N_thin
        probbase.gibbs <- array(probbase.gibbs, dim = c(N_thin, S.j, C.j))
        counter <- counter + S.j * C.j * N_thin
    }else{
        levels.gibbs <- fit[counter:(counter + N_level.j * N_thin - 1)]
        levels.gibbs <- matrix(levels.gibbs, nrow = N_thin, ncol = N_level.j, 					   byrow = TRUE)
        counter <- counter + N_level.j * N_thin
    }

    # find last time configurations
    mu_last <- matrix(fit[counter:(counter + N_sub.j * C.j - 1)], 
                     nrow = N_sub.j, ncol = C.j, byrow = TRUE)
    counter <- counter + N_sub.j * C.j
    sigma2_last <-fit[counter:(counter + N_sub.j - 1)]
    counter <- counter + N_sub.j  
    theta_last <-  matrix(fit[counter:(counter + N_sub.j * C.j - 1)], 
                         nrow = N_sub.j, ncol = C.j, byrow = TRUE)
    counter <- counter + N_sub.j * C.j     	
    

    out <- list(csmf.sub = csmf.sub, p.hat = p.hat, p.indiv = p.indiv, 
                probbase.gibbs = probbase.gibbs, 
                levels.gibbs = levels.gibbs, 
                mu.last = mu_last, sigma2.last = sigma2_last, 
                theta.last = theta_last)

    return(out)
}

# Parse_Add_Result <- function(N_sub.j, C.j, S.j, N_level.j, pool.j, fit, fit.add){
# ###########################################################
# # function to add Java output to previous results and parse into correct place 
# # @param:
# #	   various java arguments
# #	   both new and old Java output
# # @values:
# #		list of variables parsed 
#     # remove the first row from the added matrix, since it is the starting point
#     fit.all <- rbind(fit, fit.add[-1, ])
#     # parse all the saved iterations
#     out_all <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit.all)
#     # parse only the added fit to find the most recent values
#     out_last <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit.add)
#     # replace the most recent values with the correct set
#     out_all$mu.last <- out_last$mu.last
#     out_all$sigma2.last <- out_last$sigma2.last
#     out_all$theta.last <- out_last$theta.last     
#     return(out_all)
# }
#############################################################################
## InSilico VA -  helper functions for overriding standard setup
##
## author: Richard Li 
## date: 09/13/2015
#############################################################################
superInterVA.table <- function(table.dev){	
	# avoid zero
	if(min(table.dev) == 0)
	table.dev[which.min(table.dev)] <- sort(table.dev, decreasing=F)[2]/10
	return(sort(table.dev, decreasing = TRUE))
}
# superscale.vec.inter <- function(table.dev, aaa, scale = NULL, scale.max = NULL){
# ###########################################################
# # function to map ordered input vector to InterVA alphabetic scale
# # Note:  to avoid 0, change the last one to 0.000001 here
# # @param:
# #	    aaa       : vector to be scaled
# # 		scale     : sum of the vector after scaling
# #       scale.max : max of the vector after scaling
# #		reverse   : whether the order should be reversed
# # @values:	
# #  
# 	dist <- superInterVA.table(table.dev)
# 	if(length(aaa) != length(dist)){
# 		warning("dimension of probbase prior not correct, use linear prior")
# 		dist <- seq(1, 1e-5, len = length(aaa))
# 	}
# 	bbb <- dist[order(aaa)]
# 	if(!is.null(scale)) return(bbb/sum(bbb) * scale)
# 	if(!is.null(scale.max)) return(bbb * scale.max / max(bbb))
# }
superchange.inter <- function(x, table.dev, order = FALSE){
	# function to translate alphebatic matrix into numeric matrix or order matrix
	# @param:
	# 	x      : alphabetic matrix 
	#	order  : whether to change the matrix into order matrix
	# @values:
	#	numeric matrix by InterVA probbase rules, or the order matrix
		a <- dim(x)[1]
		b <- dim(x)[2]
		if(is.null(a)){
			y <- rep(0,length(x))
		}else{
			y <- matrix(0, a, b)
		}  	
		inter.table <- superInterVA.table(table.dev)
		y[x == "I"] <- inter.table[1]
	    y[x == "A+"] <- inter.table[2]
	    y[x == "A"] <- inter.table[3]
	    y[x == "A-"] <- inter.table[4]
	    y[x == "B+"] <- inter.table[5]
	    y[x == "B"] <- inter.table[6]
	    y[x == "B-"] <- inter.table[7]
	    y[x == "B -"] <- inter.table[7]
	    y[x == "C+"] <- inter.table[8]
	    y[x == "C"] <- inter.table[9]
	    y[x == "C-"] <- inter.table[10]
	    y[x == "D+"] <- inter.table[11]
	    y[x == "D"] <- inter.table[12]
	    y[x == "D-"] <- inter.table[13]
	    y[x == "E"] <- inter.table[14]
	    y[x == "N"] <- inter.table[15]
	    y[x == ""] <- inter.table[15]

	    if(order){
	    	for(i in 1: length(inter.table)){
				y[y == inter.table[i] ] <- i
			}
	    }
	    if(!is.null(a)){
			y <- matrix(y, a, b)
		}  	
	    return(y) 
}

##---------------------------------------------------------------------------------##
#############################################################################
#############################################################################
##				Helper functions all loaded                                ##
#############################################################################
#############################################################################
##---------------------------------------------------------------------------------##
	time0 <- Sys.time()
	# method <- tolower(method)
	# restrict to only normal model for now
	method <- "normal"
	alpha.scale = NULL

	## check model-specific parameters are provided
	if(method == "dirichlet"){
		if(is.null(alpha.scale)){stop("No alpha provided for Dirichlet prior")}
	}else if(method == "normal"){
		if(is.null(jump.scale)){stop("No jump scale provided for Normal prior")}
	}else{
		stop("invalid method")
	}
	if(is.null(length.sim) || is.null(thin) || is.null(burnin)){
		stop("Length of chain/thinning/burn-in not specified")
	}
	if(keepProbbase.level && Probbase_by_symp.dev){
		stop("keepProbbase.level and Probbase_by_symp.dev cannot be set to TRUE simultaneously.")
	}
	if(!is.null(nlevel.dev)){
		nlevel <- nlevel.dev
	}else{
		nlevel <- 15
	}
##---------------------------------------------------------------------------------##
## initialize key data dependencies
##	
	data("probbase", envir = environment())
	probbase<- get("probbase", envir  = environment())
	data("causetext", envir = environment())
	causetext<- get("causetext", envir  = environment())
	
	if(!dev.customization){
		# get interVA probbase
	  	prob.orig <- probbase[2:246,17:76]
	  	
	  	# get subpopulation if it's a columnname
	  	if(class(subpop) == "list" || length(subpop) == 1){
	  		col.index <- match(subpop, colnames(data))
	  		if(length(which(is.na(col.index))) > 0){
	  			stop("error: invalid sub-population name specification")
	  		}
	  		if(length(col.index) == 1){
	  			subpop <- data[, col.index]
	  		}else{
	  			subpop <- data[, col.index]
	  			subpop <- apply(subpop, 1, function(x){paste(x, collapse = ' ')})
	  		}
	  	}
	  	if(!is.null(subpop)){
	  		subpop <- as.character(subpop)
	  	} 
	  	if(length(unique(subpop)) == 1){
	  			subpop <- NULL
	  			warning("Only one level in subpopulation, running the dataset as one population")
	  	}

	  	if(dim(data)[2] != dim(probbase)[1] ){
	  		correct_names <- probbase[2:246, 2]
	  		exist <- correct_names %in% colnames(data)
	  		if(length(which(exist == FALSE)) > 0){
		        stop(paste("error: invalid data input format. Symptom(s) not found:", correct_names[!exist]))
	  		}else{
	  			data <- data[, c(1, match(correct_names, colnames(data)))]
	  		}
	    }

	    ## check the column names and give warning
	    data("RandomVA1", envir = environment())
	    RandomVA1 <- get("RandomVA1", envir  = environment())
	    valabels <- colnames(RandomVA1$data)
	    vacauses <- causetext[4:63,2]
	    external.causes = seq(41, 51)
	    external.symps = seq(211, 222)
	    
	    count.changelabel = 0
	    for(i in 1:246){
	        if(tolower(colnames(data)[i]) != tolower(valabels)[i]){
	            warning(paste("Input columne '", colnames(data)[i], "' does not match InterVA standard: '", 
	                    valabels[i], "'", sep = ""),
	                    call. = FALSE, immediate. = TRUE)
	            count.changelabel = count.changelabel + 1
	        }         
	    }
	    if(count.changelabel > 0){
	        warning(paste(count.changelabel, "column names changed in input. \n If the change in undesirable, please change in the input to match standard InterVA4 input format.\n"), call. = FALSE, immediate. = TRUE)
	        colnames(data) <- valabels
	    }
	##---------------------------------------------------------------------------------##
	  	if(!is.null(cond.prob.touse)){
	  		prob.orig <- cond.prob.touse
	  	}
	 
		#############################################################
		## remove bad data happens before taking into missing
		## (bad data refers to data without age/sex or has no real symptoms)
		tmp <- removeBad(data, isNumeric, subpop)
	  	data <- tmp[[1]]
	  	subpop <- tmp[[2]]
  	  	if(is.null(subpop)){
	  		subpop_order_list <- NULL
	  	}else{
		  	subpop_order_list <- sort(unique(subpop))
	  	}

  	}else{
  		prob.orig <- probbase.dev
		valabels <- colnames(data)
	    vacauses <- gstable.dev
  	}
  	#############################################################
  	## remove external causes
  	if(external.sep){
  		externals <- removeExt(data,prob.orig, isNumeric, subpop, subpop_order_list, external.causes, external.symps)
  		data <- externals$data
  		subpop <- externals$subpop
  		prob.orig <- externals$prob.orig
  	}
##---------------------------------------------------------------------------------##
   	## check the missing list
   	##	 note: this step is after removing bad data and before data-checking
   	## Todo: output the warnings to the warning file
   	missing.all <- NULL
   	if(TRUE){
   		## add the all missing items that are not in missing list
   		for(i in 1:(dim(data)[2]-1)){
  			if(length(which(data[,i+1] == ".")) >= 1 * dim(data)[1] && 
  			   !((i+1) %in% missing.all)){
  				missing.all <- sort(c(missing.all, i))
  			}
  		}
  		warning(paste(length(missing.all), "symptom missing completely and added to missing list", 
  			"\nList of missing symptoms: \n", 
  			paste( probbase[missing.all + 1, 2], collapse = ", ")), 
  		call. = FALSE, immediate. = TRUE)
   }
	## remove all missing symptoms from both data and probbase
	if(!is.null(missing.all)){	
		data <- data[, -(missing.all + 1)]
		prob.orig <- prob.orig[-missing.all, ]	
	}	
##---------------------------------------------------------------------------------##
	if(!dev.customization){
	  	## convert original probbase into order matrix
	  	prob.order <- change.inter(prob.orig, order = TRUE)
	  	## translate original probbase into InterVA interpreted values
	  	if(!is.numeric(prob.orig)){
	  		cond.prob.true <- change.inter(prob.orig, order = FALSE)
	 	}else{
	 		cond.prob.true <- prob.orig
	 	}
	}else{
	 	if(!is.null(prob.order.dev)){
  			prob.order <- prob.order.dev
  			cond.prob.true <- prob.orig
  		}else{
	  		prob.order <- superchange.inter(prob.orig, table.dev, order = TRUE)
	  		cond.prob.true <- superchange.inter(prob.orig, table.dev, order = FALSE)
  		}
	 } 	
##---------------------------------------------------------------------------------##
  	# get data dimensions
  	## number of data
  	N <- dim(data)[1]
	if(N <= 1) stop("Not enough sample")
	## number of symptoms (need to delete ID column)
	S <- dim(data)[2] - 1
	if(S != dim(cond.prob.true)[1]) stop("Length of symptoms is not right")
	## number of causes
	C <- dim(cond.prob.true)[2]

	## update 09/19/2015
	## external causes have been removed before here...
	if(external.sep){
		vacauses.current <- vacauses[-external.causes]
	}else{
		vacauses.current <- vacauses
	}
  if(!is.null(phy.debias)){
    phy.code <- phy.debias$code.debias
  }
	if(!is.null(phy.code)){
		#TODO: make assignment sum up to 1, and first column is unknown

		# also remove external causes too if needed
		if(external.sep){
			external.match <- match(phy.cat[, 1], vacauses[external.causes])
			phy.cat <- phy.cat[-which(!is.na(external.match)), ]
			phy.code <- phy.code[, -which(colnames(phy.code) == phy.external)]
		}
		# check every cause in phy.cat is an actual cause
		testmatch <- match(phy.cat[, 1], vacauses.current)
		if(length(which(is.na(testmatch))) > 0){
			stop("Cause name incorrect in physician category matching")
		}
		if(length(unique(phy.cat[,1])) != length(phy.cat[,1])){
			stop("Repeated cause names in physician category matching")
		}

		if(!(phy.unknown %in% colnames(phy.code))){
			stop("Cannot find Unknown category in phy.code")
		}		
		# causes in the physician coded categories
		cause.phy <- unique(phy.cat[, 2])
		if(phy.unknown %in% cause.phy){
			stop("Unknown cause exist in phy.cat matrix! Please remove it")
		}
		# Added varible for java:
		#   number of categories provided by physician + unknown
		C.phy <- length(cause.phy) + 1
		# when the two dimensions disagree
		if(dim(phy.code)[2] != C.phy + 1){
			# exist cause not appeared in phy.code, fine, just add empty column
			morecause <- which(!(cause.phy %in% colnames(phy.code)))
			if(length(morecause) > 0){
				empty <- matrix(0, dim(phy.code)[1], length(morecause))
				colnames(empty) <- cause.phy[morecause]
				phy.code <- cbind(phy.code, empty)
			}
			# exist cause not appeared in phy.cat, error
			if(dim(phy.code)[2] != C.phy + 1){
				stop("List of physician coded causes in phy.code and phy.cat do not match")
			}
		}
		# make unknown cause the first column
		matchunknown <- which(colnames(phy.code) == phy.unknown)
		restcolumns <- match(cause.phy, colnames(phy.code))
		phy.code <- phy.code[, c(1, matchunknown, restcolumns)]

		# Added varible for java:
		#   ordered set of broader cause
		vacauses.broader <- phy.cat[match(vacauses.current, phy.cat[,1]), 2]
		vacauses.broader <- match(vacauses.broader, cause.phy)
		# match ID
		matchid <- match(phy.code[, 1], data[, 1])
		# Added varible for java:
		#    assigned cause
		assignment <- matrix(0, N, C.phy)	
		# first column is unknown	
		assignment[, 1] <- 1
		# remove unmatched physician coding
		if(length(which(is.na(matchid))) > 0){
			phy.code <- phy.code[-which(is.na(matchid)), ]
			matchid <- matchid[-which(is.na(matchid))]
		}
		assignment[matchid, ] <- phy.code[, -1]
		
		#normalize assignment
		for(index in 1:dim(assignment)[1]){
			assignment[index, ] <- assignment[index, ] / sum(assignment[index, ])
		}
		
		if(external.sep){
			cat(paste(length(matchid), 
				"deaths found known physician coding after removing deaths from external causes.\n"))			
		}else{
			cat(paste(length(matchid), 
 				"deaths found known physician coding.\n"))
		}
	}else{
		# if no physician coding, everything is unknown
		C.phy <- 1
		assignment <- matrix(0, N, C.phy)	
		assignment[, 1] <- 1
		vacauses.broader <- 1:length(vacauses.current)
	}
	
##---------------------------------------------------------------------------------##
	## Specify the prior for truncated beta distribution
	prior.b.cond = trunc(1.5 * N)			
	#prior.b.cond = N

	if(keepProbbase.level){
		# if update only table of interpretation, need stronger prior
		levelcount <- table(cond.prob.true)
		# levels.strength should be multiplied here to avoid truncating to 0
		prior.b.cond <- trunc(prior.b.cond * median(levelcount) * levels.strength) 
	}else{
		prior.b.cond <- trunc(prior.b.cond * levels.strength)
	}

	if(is.null(levels.prior)){
			levels.prior <- scale.vec.inter(seq(1,nlevel), 
					scale.max = prior.b.cond * 0.99)		
	}
##---------------------------------------------------------------------------------##
	## get sub-population information
	if(!is.null(subpop)){
		if(length(subpop) != N) {
			stop("Sub-population size not match")
		}
		subpop.numeric <- rep(0, length(subpop))
		sublist <- vector("list", length(subpop_order_list))
		subbelong <- rep(0, N)
		for(i in 1:length(subpop_order_list)){
			sublist[[i]] <- which(subpop == subpop_order_list[i])
			subbelong[which(subpop == subpop_order_list[i])] <- i
			subpop.numeric[sublist[[i]]] <- i-1
		}
		names(sublist)<- subpop_order_list
		N.sub <- length(sublist)
	}
##---------------------------------------------------------------------------------##
  	## initiate indicator matrix	
	indic <- matrix(0, N, S)
	id <- data[, 1]
	## convert "Y" to 1 in indicator matrix
	## if numeric data, just use data metrix without ID
	if(isNumeric){
		indic <- data[, -1]
	}else{		
		for(i in 1:S){
			temp <- toupper(data[,i+1])
			temp.ind <- which(temp == "Y")
			indic[temp.ind, i] <- 1
		}
	}
	## if data contains missing ".", translate to -1
	if(isNumeric){
		contains.missing <- (length(which(data == "-1"))> 0)
	}else{
		contains.missing <- (length(which(data == ".")) > 0)
	}
	# if there are missing and data not numeric, initiate missing
	if(contains.missing && !(isNumeric)){
		for(i in 1:S){
			temp <- data[,i+1]
			temp.ind <- which(temp == ".")
			indic[temp.ind, i] <- -1
		}
	}
	## if data contains subpopulation, 
	## for complete missing within subpopulation,
	##	   change missing from -1 to -2 to separate.
	if(!is.null(subpop)){
		for(j in 1:N.sub){
			for(i in 1:S){
				temp <- indic[which(subbelong == j), i]
				if(length(which(temp == -1)) == length(temp)){
					indic[which(subbelong == j), i] <- -2	
				}
			}			
		}
	}
##---------------------------------------------------------------------------------##
	if(datacheck){
		cat("Performing data consistency check...\n")
		indic <- datacheck.interVA(id, indic, missing.all, external.sep, warning.write)
		cat("Data check finished.\n")
	}
	indic.w.missing <- indic
##---------------------------------------------------------------------------------##
## parameter initialization
   # csmf.prior <- rep(1/C, C)

	Sys_Prior <- as.numeric(change.inter(probbase[1,17:76], order = FALSE))
	# Number of indicators + 13 description variables. A_group:14-16;B_group:17:76;D_group:77:81
	D <- length(Sys_Prior)
	csmf.prior <- Sys_Prior/sum(Sys_Prior)
##---------------------------------------------------------------------------------##
##  initialize prior and adjust for external causes
	if(method == "dirichlet"){
		# determine alpha for CSMF
  		alpha <- csmf.prior * C * alpha.scale
  	}else if(method == "normal"){
  		if(!is.null(csmf.prior)){
	  		Z <- csmf.prior[1]/exp(1)
	  		mu <- log(csmf.prior / Z)
	  	}else{
	  		mu <- rep(1, C)
	  	}
	  	sigma2 <- 1
  	}
	if(external.sep){
		csmf.prior <- csmf.prior[-external.causes]
		if(method == "normal" && length(mu) != C){ 
			mu <- mu[-external.causes]
		}else if(method == "dirichlet" && length(alpha) != C){
			alpha <- alpha[-external.causes]
		}
	}
##---------------------------------------------------------------------------------##
## initialize probbase and start java
	
	if(!useProbbase){
		cond.prob <- cond.initiate(prob.order, expIni = TRUE, Inter.ini = TRUE,
						  min = trunc.min, max = trunc.max)
    }else{
    	cond.prob <- change.inter(prob.orig)
    }

    # library(rJava)
	if(is.null(java_option)) java_option = "-Xmx1g"
	options( java.parameters = java_option )

	obj <- .jnew("sampler/InsilicoSampler2")

    N.j <- as.integer(N)
    S.j <- as.integer(S)
    C.j <- as.integer(C)
    probbase.j <- .jarray(as.matrix(cond.prob), dispatch=TRUE)
    
     if(dev.customization && !(useProbbase)){
    	# get new dist
    	dist <- superInterVA.table(table.dev)
    	# check existence of probbase levels
    	exist <- seq(1:length(dist)) %in% unique(as.vector(prob.order))
    	
    	# update order matrix
    	prob.order.new <- prob.order
    	for(i in 1:length(dist)){
    		if(!exist[i]){
    			# I is 1, N is 15
    			# if a level is missing, all level after should minus 1
    			prob.order.new[which(prob.order > i)] <- prob.order.new[which(prob.order > i)] - 1
    		}
    	}
    	prob.order <- prob.order.new
    	probbase_order.j <- .jarray(as.matrix(prob.order), dispatch = TRUE)
    	# update level vector
    	dist <- dist[exist]
    	levels.prior <- levels.prior[exist]
    	N_level.j <- as.integer(sum(exist))
    }else{
    	probbase_order.j <- .jarray(as.matrix(prob.order), dispatch = TRUE)
    	N_level.j <- as.integer(nlevel)
    	dist <- InterVA.table(0)
	}
   
    level_values.j <- .jarray(dist, dispatch = TRUE)
    prior_a.j <- .jarray(levels.prior , dispatch = TRUE)
    prior_b.j <- prior.b.cond
    jumprange.j <- jump.scale
    trunc_min.j <- trunc.min
    trunc_max.j <- trunc.max
    indic.j <- .jarray(as.matrix(indic), dispatch=TRUE)
    contains_missing.j <- as.integer(contains.missing)
    # update Oct 1, 2015: pool variable redefined:
    # 					  0 - pool to table; 1 - by cause; 2 - by symptom
    pool.j <- as.integer(!keepProbbase.level) + as.integer(Probbase_by_symp.dev)
    seed.j <- as.integer(seed)
    N_gibbs.j <- as.integer(length.sim)
    burn.j <- as.integer(burnin)
    thin.j <- as.integer(thin)
    mu.j <- .jarray(mu, dispatch = TRUE)
    sigma2.j <- sigma2 
    isUnix <-  .Platform$OS.type == "unix"

    assignment.j <- .jarray(as.matrix(assignment), dispatch = TRUE)
    C.phy.j <- as.integer(C.phy)
    vacauses.broader.j <- .jarray(vacauses.broader+0.0, dispatch = TRUE)

    if(is.null(subpop)){
		N_sub.j <- as.integer(1)
		subpop.j <- .jarray(as.integer(rep(0, N)), dispatch = TRUE)
	}else{
	    N_sub.j <- as.integer(N.sub)
        subpop.j <- .jarray(as.integer(subpop.numeric), dispatch = TRUE)
	}

    isAdded <- FALSE
    mu.last.j <- .jarray(matrix(0, N_sub.j, C), dispatch = TRUE)
    sigma2.last.j <- .jarray(rep(0, N_sub.j), dispatch = TRUE)
	theta.last.j <- .jarray(matrix(0, N_sub.j, C), dispatch = TRUE)

    ins  <- .jcall(obj, "[D", "Fit", 
		N.j, S.j, C.j, N_sub.j, N_level.j, 
		probbase.j, probbase_order.j, level_values.j, 
		prior_a.j, prior_b.j, jumprange.j, trunc_min.j, trunc_max.j, 
		indic.j, subpop.j, contains_missing.j, pool.j, 
		seed.j, N_gibbs.j, burn.j, thin.j, 
		mu.j, sigma2.j, isUnix, useProbbase, 
		isAdded, mu.last.j, sigma2.last.j, theta.last.j, 
		C.phy.j, vacauses.broader.j, assignment.j) 
    # one dimensional array is straightforward
    fit <- ins
    # fit <-  .jeval(ins, .jevalArray))
##---------------------------------------------------------------------------------##	
##
## save data from java output to proper format
##    
    results <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit)
	# check convergence
    conv <- tryCatch({
					   if(!is.null(results$csmf.sub)){
							csmf.diag(results$csmf.sub, conv.csmf, 
								              test = "heidel", verbose = FALSE) 	
				    	}else{
				    		csmf.diag(results$p.hat, conv.csmf, 
								              test = "heidel", verbose = FALSE) 
				    	}
				}, error = function(condition) {
				    print("error checking diagnostics")
				    FALSE
				})
    # check convergence and if the chain needs to run longer
    if(auto.length){	
    	# if not converge, run again, max number of runs = 3
    	add <- 1
    	while(!conv && add < 3){
    		mu.last.j <- .jarray(as.matrix(results$mu.last), dispatch = TRUE)
    		sigma2.last.j <- .jarray(results$sigma2.last, dispatch = TRUE)
    		theta.last.j <- .jarray(as.matrix(results$theta.last), dispatch = TRUE)
    		# same length as previous chain if added the first time
    		# double the length	if the second time
    		length.sim <- length.sim * 2 
    		burnin <- length.sim / 2
    		N_gibbs.j <- as.integer(trunc(N_gibbs.j * (2^(add-1))))
			burn.j <- as.integer(0)

			cat(paste("Not all causes with CSMF >", conv.csmf, "are convergent.\n"))
    		cat(paste("Increase chain length with another", N_gibbs.j, "iterations\n"))
    		obj <- .jnew("sampler/InsilicoSampler2")
    		ins  <- .jcall(obj, "[D", "Fit", 
						N.j, S.j, C.j, N_sub.j, N_level.j, 
						probbase.j, probbase_order.j, level_values.j, 
						prior_a.j, prior_b.j, jumprange.j, trunc_min.j, trunc_max.j, 
						indic.j, subpop.j, contains_missing.j, pool.j, 
						seed.j, N_gibbs.j, burn.j, thin.j, 
						mu.j, sigma2.j, isUnix, useProbbase, 
						TRUE, mu.last.j, sigma2.last.j, theta.last.j, 
						C.phy.j, vacauses.broader.j, assignment.j)
    		# one dimensional array is straightforward
    		fit.add <- ins
    		# fit.add <-  t(sapply(ins, .jevalArray))

    		# remove the first row of fit)
    		results <- ParseResult(N_sub.j, C.j, S.j, N_level.j, pool.j, fit.add)
    		add = add + 1
    		# check convergence
	    	if(!is.null(results$csmf.sub)){
				conv <- csmf.diag(results$csmf.sub, conv.csmf, 
					              test = "heidel", verbose = FALSE) 	
	    	}else{
	    		conv <- csmf.diag(results$p.hat, conv.csmf, 
					              test = "heidel", verbose = FALSE) 
	    	}
    	}
    }
    ## if still not convergent 
    if(!conv){
    	cat(paste("Not all causes with CSMF >", conv.csmf, "are convergent.\n",
    			  "Please check using csmf.diag() for more information.\n"))
    }
    csmf.sub  <- results$csmf.sub 
    p.hat  <- results$p.hat
    p.indiv  <- results$p.indiv 
    probbase.gibbs  <- results$probbase.gibbs 
    levels.gibbs  <- results$levels.gibbs
##---------------------------------------------------------------------------------##
## To make output consistent for further analysis,
## 		add back external results
	##
	## if separated external causes
	##
    if(external.sep){
    	# get starting and ending index
    	ext1 <- external.causes[1]
    	ext2 <- external.causes[2]
    	##
    	## if with subgroup
    	##
    	if(!is.null(subpop)){
    		# set p.hat to NULL, and set up csmf.sub.all as a list of p.hat 
    		p.hat <- NULL
    		csmf.sub.all <- vector("list", N_sub.j)
    		names(csmf.sub.all) <- subpop_order_list
    		# iterate over all subpopulation
    		for(j in 1:length(csmf.sub)){
    			# initialize the csmf matrix 
    			csmf.sub.all[[j]] <- matrix(0, 
    				dim(csmf.sub[[j]])[1], 
    				C.j + length(external.causes))
    			# rescale the non-external CSMF once the external causes are added 
    			rescale <- length(sublist[[j]]) / (length(sublist[[j]]) + length(which(externals$ext.sub == subpop_order_list[j])))
    			temp <- csmf.sub[[j]] * rescale
    					
    			# combine the rescaled non-external CSMF with the external CSMF	
    			csmf.sub.all[[j]]  <- cbind(temp[, 1:(ext1 - 1)], 
				    						matrix(externals$ext.csmf[[j]], 
				    						  	   dim(temp)[1],
				    						  	   length(external.causes), 
				    						  	   byrow = TRUE), 
				    						temp[, ext1:C.j])	
    		}
    		csmf.sub <- csmf.sub.all
    	##
    	## if no subgroup
    	##
    	}else{
    		p.hat <- p.hat * N/(N + length(externals$ext.id))
    		temp <- p.hat[, ext1:C.j]
    		extra <- matrix(externals$ext.csmf, dim(p.hat)[1], length(external.causes), byrow = TRUE)
    		p.hat <- cbind(p.hat[, 1:(ext1 - 1)], extra, temp)
    	}

    	p.indiv <- cbind(p.indiv[, 1:(ext1 - 1)], 
    						  matrix(0, dim(p.indiv)[1], length(external.causes)), 
    						  p.indiv[, (ext1:C.j)])
    	
    	p.indiv.ext <- matrix(0, nrow = length(externals$ext.id), ncol = C.j + length(external.causes) )
    	for(i in 1:length(externals$ext.id)){p.indiv.ext[i, externals$ext.cod[i]] <- 1}
    	p.indiv <- rbind(p.indiv, p.indiv.ext) 
    	id <- c(id, externals$ext.id)
    	subpop <- c(subpop, externals$ext.sub)
   }
##---------------------------------------------------------------------------------##    	
## add column names to outcome
if(pool.j != 0){
	if(external.sep){
		# remove column "ID"
		valabels <- valabels[-1]
		# remove external
		valabels <- valabels[-external.symps]
		vacauses.ext <- vacauses[-external.causes]
	    # remove missing
		valabels <- valabels[-missing.all]
	}
	dimnames(probbase.gibbs)[[2]] <- valabels
	dimnames(probbase.gibbs)[[3]] <- vacauses.ext		
}else{
	if(dev.customization && is.null(nlevel.dev)){
		colnames(levels.gibbs) <- c("I", "A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D+", "D", "D-", "E", "N")[exist]
	}else{
		colnames(levels.gibbs) <- c("I", "A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D+", "D", "D-", "E", "N")
	}
	probbase.gibbs <- levels.gibbs
} 
if(!is.null(subpop)){
	  for(j in 1:length(csmf.sub)){colnames(csmf.sub[[j]]) <- vacauses}
	  p.hat <- csmf.sub
}else{
	  colnames(p.hat) <- vacauses
}
colnames(p.indiv) <- vacauses
rownames(p.indiv) <- id
##---------------------------------------------------------------------------------##    	
if(useProbbase){
    	probbase.gibbs <- NULL
}
out <- list(
		id = id,
	    indiv.prob = p.indiv, 
		csmf = p.hat,
		conditional.probs = probbase.gibbs,
		missing.symptoms = missing.all,
		external = external.sep, 
		
		useProbbase = useProbbase, 
		keepProbbase.level = keepProbbase.level, 
		datacheck = datacheck,
		length.sim = length.sim, 
		thin = thin, 
		burnin = burnin, 
	
		jump.scale = jump.scale, 
		levels.prior = levels.prior, 
		levels.strength = levels.strength, 
		trunc.min = trunc.min, 
		trunc.max = trunc.max, 
		subpop = subpop)
class(out) <- "insilico"		
return(out)  	
} 
