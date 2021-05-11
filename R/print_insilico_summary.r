#' Print method for summarizing InSilicoVA Model Fits
#' 
#' This function is the print method for class \code{insilico_summary}.
#' 
#' 
#' @param x \code{insilico_summary} object.
#' @param ... not used
#' 
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{summary.insilico}} 
#' @references 
#' Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C. Crampin,
#' Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death assignment
#' using verbal autopsies, \emph{Journal of the American Statistical
#' Association} (2016), 111(515):1036-1049.
#' @examples
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
#'               Nsim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE)
#' summary(fit1)
#' summary(fit1, top = 10)
#'
#' # save individual COD distributions to files
#' summary(fit1, file = "results.csv")
#' }
#' @export 
print.insilico_summary <- function(x, ...) {
	# print single death summary
	if(!is.null(x$indiv.prob)){
		cat(paste0("InSilicoVA fitted top ", x$top, " causes for death ID: ", x$id, "\n"))		
		cat(paste0("Credible intervals shown: ", round(x$indiv.CI * 100), "%\n"))
		print(round(x$indiv.prob, 4), digits = 4)

	# print population summary
	}else{
		cat("InSilicoVA Call: \n")
		cat(paste(length(x$id.all), "death processed\n"))
		cat(paste(x$Nsim, "iterations performed, with first", 
				  x$burnin, "iterations discarded\n",
				  trunc((x$Nsim - x$burnin)/x$thin), "iterations saved after thinning\n"))
		if(!x$updateCondProb){
				cat("Fitted with fixed conditional probability matrix\n")
			}else if(x$keepProbbase.level){
				cat("Fitted with re-estimated conditional probability level table\n")	
			}else{
				cat("Fitted with re-estimating conditional probability matrix\n")	
			}   
		if(x$datacheck){
			interva.version <- ifelse(x$data.type == "WHO2012", "InterVA4", "InterVA5")
			cat(paste0("Data consistency check performed as in ", interva.version,  "\n"))
		}

		if(!is.null(x$subpop_counts)){
			cat("Sub population frequencies:")
			print(x$subpop_counts)
		}	
		
		if(!methods::is(x$csmf.ordered, "list")){
			top <- min(x$showTop, dim(x$csmf.ordered)[1])
			cat("\n")
			cat(paste("Top", top,  "CSMFs:\n"))
			csmf.out.ordered <- x$csmf.ordered
		    out <- as.matrix(csmf.out.ordered[1:top, ])
				out[, 1] <- round(out[, 1], 4)
				out[, 2] <- round(out[, 2], 4)
				out[, 3] <- round(out[, 3], 4)
				out[, 4] <- round(out[, 4], 4)
				out[, 5] <- round(out[, 5], 4)
			print(out)		
		}else{
			for(i in 1:length(x$csmf.ordered)){
				top <- min(x$showTop, dim(x$csmf.ordered[[i]])[1])
				cat("\n")
				cat(paste(names(x$csmf.ordered)[i], "- Top", top,  "CSMFs:\n"))
				csmf.out.ordered <- x$csmf.ordered[[i]]
			    out <- as.matrix(csmf.out.ordered[1:top, ])
					out[, 1] <- round(out[, 1], 4)
					out[, 2] <- round(out[, 2], 4)
					out[, 3] <- round(out[, 3], 4)
					out[, 4] <- round(out[, 4], 4)
					out[, 5] <- round(out[, 5], 4)
				print(out)		
			}
		}
	}
	
}
	

