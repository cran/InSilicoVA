
print.insilico_summary <- function(x, ...) {
	cat("InSilicoVA Call: \n")
	cat(paste(length(x$id), "death processed\n"))
	cat(paste(x$length.sim, "iterations performed, with first", 
			  x$burnin, "iterations discarded\n",
			  trunc((x$length.sim - x$burnin)/x$thin), "iterations saved after thinning\n"))
	if(x$useProbbase){
			cat("Fitted with InterVA4 conditional probability matrix\n")
		}else if(x$keepProbbase.level){
			cat("Fitted with re-estimated InterVA4 conditional probability level table\n")	
		}else{
			cat("Fitted with re-estimating InterVA4 conditional probability matrix\n")	
		}   
	if(x$datacheck){
		cat("Data consistency check performed as in InterVA4 \n")
	}
	if(!is.null(x$HIV)){
		cat(paste("HIV level:     ", x$HIV, "\n"))
		cat(paste("Malaria level: ", x$Malaria, "\n"))
	}

	if(!is.null(x$subpop_counts)){
		cat("Sub population frequencies:")
		print(x$subpop_counts)
	}	
	
	if(class(x$csmf.ordered) != "list"){
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
	

