insilico_v0 <- function(data, isNumeric = FALSE,useProbbase = FALSE, keepProbbase.level = TRUE,  cond.prob.touse = NULL,datacheck = TRUE, warning.write = FALSE, external.sep = TRUE, length.sim = 4000, thin = 10, burnin = 2000, auto.length = TRUE, conv.csmf = 0.02, jump.scale = 0.1, levels.prior = NULL, levels.strength = 1, trunc.min = 0.0001, trunc.max = 0.9999, subpop = NULL, java_option = "-Xmx1g", seed = 1){ 
	
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
  
removeExt <- function(data, prob.orig, is.Numeric, subpop, external.causes, external.symps){
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
		ext.csmf <- vector("list", length(unique(subpop)))
		for(i in 1:length(ext.csmf)){
			ext.csmf[[i]] <- rep(0, length(extCauses))
			ext.cod.temp <- ext.cod[which(ext.sub == unique(subpop)[i])]
			if(!is.null(ext.cod.temp)){
				for(j in 1:length(extCauses)){
					ext.csmf[[i]][j] <- length(which(ext.cod.temp == extCauses[j]))
				}
				ext.csmf[[i]] <- ext.csmf[[i]]/length(which(subpop == unique(subpop)[i]))		
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
    if(pool.j == 0){
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

##---------------------------------------------------------------------------------##
## initialize key data dependencies
##	
	data("probbase", envir = environment())
	probbase<- get("probbase", envir  = environment())
	data("causetext", envir = environment())
	causetext<- get("causetext", envir  = environment())
	
	# get interVA probbase
  	prob.orig <- probbase[2:246,17:76] 

  	if(dim(data)[2] != dim(probbase)[1] ){
        stop("error: invalid data input format. Number of values incorrect")
    }

    ## check the column names and give warning
    data("SampleInput_insilico", envir = environment())
    SampleInput_insilico <- get("SampleInput_insilico", envir  = environment())
    valabels <- colnames(SampleInput_insilico$data)
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
 #  	if(useInterVA == 2){
	# 		va <- InterVA(as.matrix(data), HIV = HIV, Malaria = Malaria, write = FALSE)
	#   		csmf.inter <- Population.summary(va$VA, min.prob = 0, noplot = TRUE)
	#   		if(length(csmf.inter) == 61) csmf.inter <- as.numeric(csmf.inter[-61])
	# }	
	#############################################################
	## remove bad data happens before taking into missing
	## (bad data refers to data without age/sex or has no real symptoms)
	tmp <- removeBad(as.matrix(data), isNumeric, subpop)
  	data <- tmp[[1]]
  	subpop <- tmp[[2]]

  	#############################################################
  	## remove external causes
  	if(external.sep){
  		externals <- removeExt(data,prob.orig, isNumeric, subpop, external.causes, external.symps)
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
  	## convert original probbase into order matrix
  	prob.order <- change.inter(prob.orig, order = TRUE)
  	## translate original probbase into InterVA interpreted values
  	if(!is.numeric(prob.orig)){
  		cond.prob.true <- change.inter(prob.orig, order = FALSE)
 	}else{
 		cond.prob.true <- prob.orig
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
			levels.prior <- scale.vec.inter(seq(1,15), 
					scale.max = prior.b.cond * 0.99)		
	}
##---------------------------------------------------------------------------------##
	## get sub-population information
	if(!is.null(subpop)){
		if(length(subpop) != N) {
			stop("Sub-population size not match")
		}
		subpop.numeric <- rep(0, length(subpop))
		sublist <- vector("list", length(unique(subpop)))
		subbelong <- rep(0, N)
		for(i in 1:length(unique(subpop))){
			sublist[[i]] <- which(subpop == unique(subpop)[i])
			subbelong[which(subpop == unique(subpop)[i])] <- i
			subpop.numeric[sublist[[i]]] <- i-1
		}
		names(sublist)<- unique(subpop)
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
			temp <- data[,i+1]
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

	obj <- .jnew("sampler/InsilicoSampler")

    N.j <- as.integer(N)
    S.j <- as.integer(S)
    C.j <- as.integer(C)
    N_level.j <- as.integer(15)
    probbase.j <- .jarray(as.matrix(cond.prob), dispatch=TRUE)
    probbase_order.j <- .jarray(as.matrix(prob.order), dispatch = TRUE)

    dist <- InterVA.table(0)
    level_values.j <- .jarray(dist, dispatch = TRUE)
    prior_a.j <- .jarray(levels.prior , dispatch = TRUE)
    prior_b.j <- prior.b.cond
    jumprange.j <- jump.scale
    trunc_min.j <- trunc.min
    trunc_max.j <- trunc.max
    indic.j <- .jarray(as.matrix(indic), dispatch=TRUE)
    contains_missing.j <- as.integer(contains.missing)
    pool.j <- as.integer(keepProbbase.level)
    seed.j <- as.integer(seed)
    N_gibbs.j <- as.integer(length.sim)
    burn.j <- as.integer(burnin)
    thin.j <- as.integer(thin)
    mu.j <- .jarray(mu, dispatch = TRUE)
    sigma2.j <- sigma2 
    isUnix <-  .Platform$OS.type == "unix"

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
		isAdded, mu.last.j, sigma2.last.j, theta.last.j) 
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
    		obj <- .jnew("sampler/InsilicoSampler")
    		ins  <- .jcall(obj, "[D", "Fit", 
						N.j, S.j, C.j, N_sub.j, N_level.j, 
						probbase.j, probbase_order.j, level_values.j, 
						prior_a.j, prior_b.j, jumprange.j, trunc_min.j, trunc_max.j, 
						indic.j, subpop.j, contains_missing.j, pool.j, 
						seed.j, N_gibbs.j, burn.j, thin.j, 
						mu.j, sigma2.j, isUnix, useProbbase, 
						TRUE, mu.last.j, sigma2.last.j, theta.last.j) 
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
    		names(csmf.sub.all) <- unique(subpop)
    		# iterate over all subpopulation
    		for(j in 1:length(csmf.sub)){
    			# initialize the csmf matrix 
    			csmf.sub.all[[j]] <- matrix(0, 
    				dim(csmf.sub[[j]])[1], 
    				C.j + length(external.causes))
    			# rescale the non-external CSMF once the external causes are added 
    			rescale <- length(sublist[[j]]) / (length(sublist[[j]]) + length(which(externals$ext.sub == unique(subpop)[j])))
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
    	# following codes used when p.indiv is recording all iterations
    	# for(i in 1:length(externals$ext.id)){p.indiv.ext[i, externals$ext.cod[i]] <- 1}
    	p.indiv <- rbind(p.indiv, p.indiv.ext) 
    	id <- c(id, externals$ext.id)
   }
##---------------------------------------------------------------------------------##    	
## add column names to outcome
if(pool.j == 0){
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
	colnames(levels.gibbs) <- c("I", "A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D+", "D", "D-", "E", "N")
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
