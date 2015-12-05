##
## function to check if the csmf chains have converged
##
## csmf     : the csmf to be checked, either vector or list (with sub-population)
## conv.csmf: the minimum mean csmf to check
## test     : type of test
## verbose  : whether to return the test statistics, or simply outcome



#' Convergence test for fitted InSilico model
#' 
#' Produce convergence test for CSMFs from fitted \code{"insilico"} objects.
#' 
#' The tests are performed using \link{heidel.diag} and \link{gelman.diag}
#' functions in \code{coda} package. The function takes either one or a list of
#' output from \code{insilico} function, or only the iteration by CSMF matrix.
#' Usually in practice, many causes with very tiny CSMF are hard to converge
#' based on standard tests, thus it is suggested to check convergence for only
#' causes with mean CSMF over certain threshold by setting proper
#' \code{conv.csmf}.
#' 
#' Note for Gelman and Rubin's test, all chains should have the same length. If
#' the chains are sampled with automatically length determination, they might
#' not be comparable by this test.
#' 
#' @param csmf It could be either fitted \code{"insilico"} object, a list of
#' fitted \code{"insilico"} object from different chains of the same length but
#' different starting values, i.e., different seed. Or it could be the matrix
#' of CSMF obtained by \code{insilico}, or the list of matrices of CSMF. All
#' CSMF could contain more than one subpopulations, but should be in the same
#' format and order.
#' @param conv.csmf The minimum mean CSMF to be checked. Default to be 0.02,
#' which means any causes with mean CSMF lower than 0.02 will not be tested.
#' @param test Type of test. Currently supporting Gelman and Rubin's test
#' (\code{test = "gelman"}) for multi-chain test, and Heidelberger and Welch's
#' test (\code{test = "heidel"}) for single-chain test.
#' @param verbose Logical indicator to return the test detail instead of one
#' logical outcome for Heidelberger and Welch's test. Default to be TRUE.
#' @param autoburnin Logical indicator of whether to omit the first half of the
#' chain as burn in. Default to be FALSE since \code{insilico} return only the
#' iterations after burnin by default.
#' @param external.sep Logical indicator for whether to separate out external
#' causes first. Default set to be TRUE. If set to TRUE, the function will not
#' check convergence for external causes, e.g., traffic accident, accidental
#' fall, suicide, etc.
#' @param \dots Arguments to be passed to \link{heidel.diag} or
#' \link{gelman.diag}
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{insilico}}, \code{\link{summary.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' 
#' Gelman, Andrew, and Donald B. Rubin. Inference from iterative simulation
#' using multiple sequences. \emph{Statistical science} (1992): 457-472.
#' 
#' Brooks, Stephen P., and Andrew Gelman. General methods for monitoring
#' convergence of iterative simulations. \emph{Journal of computational and
#' graphical statistics} 7.4 (1998): 434-455.
#' 
#' Heidelberger, Philip, and Peter D. Welch. A spectral method for confidence
#' interval generation and run length control in simulations.
#' \emph{Communications of the ACM} 24.4 (1981): 233-245.
#' 
#' Heidelberger, Philip, and Peter D. Welch. Simulation run length control in
#' the presence of an initial transient. \emph{Operations Research} 31.6
#' (1983): 1109-1144.
#' 
#' Schruben, Lee W. Detecting initialization bias in simulation output.
#' \emph{Operations Research} 30.3 (1982): 569-590.
#' @keywords InSilicoVA
#' @examples
#' 
#' # load sample data together with sub-population list
#' data(RandomVA1)
#' \donttest{
#' # extract INterVA style input data
#' data <- RandomVA1$data
#' # extract sub-population information. 
#' # The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
#' subpop <- RandomVA1$subpop
#' 
#' # run without sub-population
#' fit1a<- insilico( data, subpop = NULL, 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1, 
#'               auto.length = FALSE)
#' fit1b<- insilico( data, subpop = NULL,  
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 2, 
#'               auto.length = FALSE)
#' fit1c<- insilico( data, subpop = NULL,  
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 3, 
#'               auto.length = FALSE)
#' # single chain check
#' csmf.diag(fit1a)
#' # equivalent way of check one chain
#' csmf.diag(fit1a$csmf)
#' 
#' # multiple chains check
#' csmf.diag(list(fit1a, fit1b, fit1c), test = "gelman")
#' # equivalent way of check one chain
#' csmf.diag(list(fit1a$csmf, fit1b$csmf, fit1c$csmf), test = "gelman")
#' 
#' 
#' # with sub-populations
#' fit2a<- insilico( data, subpop = subpop,  
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1, 
#'               auto.length = FALSE)
#' fit2b<- insilico( data, subpop = subpop,  
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 2, 
#'               auto.length = FALSE)
#' fit2c<- insilico( data, subpop = subpop,   
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 3, 
#'               auto.length = FALSE)
#' 
#' # single chain check
#' csmf.diag(fit2a)
#' # equivalent way of check one chain
#' csmf.diag(fit2a$csmf)
#' 
#' # multiple chains check
#' csmf.diag(list(fit2a, fit2b, fit2c), test = "gelman")
#' # equivalent way of check one chain
#' csmf.diag(list(fit2a$csmf, fit2b$csmf, fit2c$csmf), test = "gelman")
#' }
#' 
#' 
#' @export csmf.diag
csmf.diag <- function(csmf, conv.csmf = 0.02, test= c("gelman", "heidel")[2], verbose = TRUE, autoburnin = FALSE, external.sep = TRUE, ...){
	
	external.causes = seq(41, 51)
	
	# check if the input is insilico object or only csmf
	if(class(csmf) == "insilico"){
		if(external.sep){
			csmf <- csmf[, -external.causes]
		}else{
			csmf <- csmf$csmf
		}
	}  
	if(class(csmf) == "list" && class(csmf[[1]]) == "insilico"){
		csmf_only <- vector("list", length(csmf))
		for(i in 1:length(csmf)){
			if(external.sep){
				csmf_only[[i]] <- csmf[[i]]$csmf[, -external.causes]
			}else{
				csmf_only[[i]] <- csmf[[i]]$csmf
			}
		}
		csmf <- csmf_only
	}
	
	# check conv.csmf value
	if(conv.csmf >= 1 || conv.csmf <= 0) {
		stop("Invalid threshold.")
	}
	if(test == "gelman" && class(csmf) != "list"){
		stop("Need more than one chain to perform Gelman and Rubin's test.")
	}
	gelman.single <- function(list, conv.csmf){
		# initialize the mean csmf for all chains
		mean <- rep(0, dim(list[[1]])[2])
		# in case chains have different lengths, count the normalizing constant
		Z <- 0
		for(i in 1:length(list)){
			# add sums
			mean <- mean + apply(list[[i]], 2, sum)
			# add counts
			Z <- Z + dim(list[[i]])[1]
			# change to mcmc object
			list[[i]] <- mcmc(list[[i]])
		}
		# sort and find which causes to report
		mean <- mean / Z
		which <- which(mean > conv.csmf)
		order <- order(mean[which], decreasing = TRUE)
		for(i in 1:length(list)){
			list[[i]] <- list[[i]][, which[order]]
		}
		# get the mcmc list
		multi <- mcmc.list(list)
		test <- gelman.diag(multi, autoburnin=FALSE)
		return(test)
	}
	# function to perform heidel test on one chain
	heidel.single <- function(one, conv.csmf){
		one <- as.matrix(one)
		# sort chain by CSMF mean values
		mean <- apply(one, 2, mean)
		one <- one[, order(mean, decreasing = TRUE)]
		# remove CSMFs below conv.csmf
		mean <- apply(one, 2, mean)
		# use drop = FALSE to filter without making one an vector
		one <- one[, which(mean > conv.csmf), drop = FALSE]
		one <- mcmc(one)
		test <- heidel.diag(one)
		return(test)
	}
	# Heidel test
	if(test == "heidel"){
		testout <- NULL
		conv <- 1
		if(class(csmf) == "list"){
			testout <- vector("list", length(csmf))
			names(testout) <- names(csmf)
			for(i in 1:length(testout)){
				testout[[i]] <- heidel.single(csmf[[i]], conv.csmf)
				conv <- conv * prod(testout[[i]][, 1]) * prod(testout[[i]][, 4])
			}
		}else{
			testout <- heidel.single(csmf, conv.csmf)
			conv <- prod(testout[, 1]) * prod(testout[, 4])
		}
	}
	# Gelman test
	if(test == "gelman"){
		testout <- NULL
		conv <- 1
		if(class(csmf[[1]]) == "list"){
			testout <- vector("list", length(csmf[[1]]))
			names(testout) <- names(csmf[[1]])
			for(i in 1:length(testout)){
				# create new list for each sub population
				csmf_sub <- vector("list", length(csmf))
				for(j in 1:length(csmf)){
					csmf_sub[[j]] <- csmf[[j]][[i]]
				}
				testout[[i]] <- gelman.single(csmf_sub, conv.csmf)
			}
		}else{
			testout <- gelman.single(csmf, conv.csmf)
		}
	}
	
	if(is.na(conv)) conv <- 0
	if(!verbose && test == "heidel"){
		return(conv)
	}else{
		return(testout)
	}
}
