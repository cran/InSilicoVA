#' plot CSMF from a "insilico" object
#' 
#' Produce a bar plot of the CSMFs for a fitted \code{"insilico"} object.
#' 
#' To-do
#' 
#' @param x fitted \code{"insilico"} object
#' @param type An indicator of the type of chart to plot. "errorbar" for line
#' plots of only the error bars on single population; "bar" for bar chart with
#' error bars on single population; "compare" for line charts on multiple
#' sub-populations.
#' @param top The number of top causes to plot. If multiple sub-populations are
#' to be plotted, it will plot the union of the top causes in all
#' sub-populations.
#' @param causelist The list of causes to plot. It could be a numeric vector
#' indicating the position of the causes in the InterVA cause list (see
#' \code{\link{causetext}}), or a vector of character string of the cause
#' names. The argument supports partial matching of the cause names. e.g.,
#' "HIV/AIDS related death" could be abbreviated into "HIV"; "Other and
#' unspecified infect dis" could be abbreviated into "Other and unspecified
#' infect".
#' @param which.sub Specification of which sub-population to plot if there are
#' multiple and \code{type} is set to "bar".
#' @param xlab Labels for the causes.
#' @param ylab Labels for the CSMF values.
#' @param title Title of the plot.
#' @param horiz Logical indicator indicating if the bars are plotted
#' horizontally.
#' @param angle Angle of rotation for the texts on x axis when \code{horiz} is
#' set to FALSE
#' @param fill The color to fill the bars when \code{type} is set to "bar".
#' @param border The color to color the borders of bars when \code{type} is set
#' to "bar".
#' @param err_width Size of the error bars.
#' @param err_size Thickness of the error bar lines.
#' @param point_size Size of the points.
#' @param bw Logical indicator for setting the theme of the plots to be black
#' and white.
#' @param \dots Not used.
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{insilico}}, \code{\link{summary.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' @keywords InSilicoVA
#' @examples
#' 
#' \dontrun{
#' # load sample data together with sub-population list
#' data(RandomVA1)
#' # extract INterVA style input data
#' data <- RandomVA1$data
#' # extract sub-population information. 
#' # The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
#' subpop <- RandomVA1$subpop
#' 
#' # run without sub-population
#' fit1<- insilico( data, subpop = NULL, 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE, 
#'				 auto.length = FALSE)
#' # default plot
#' plot(fit1)
#' 
#' # customized line plot
#' plot(fit1, top = 15, horiz = FALSE, fill = "gold", 
#' 		   bw = TRUE, title = "Top 15 CSMFs", angle = 70, 
#'        err_width = .2, err_size = .6, point_size = 2)
#' 
#' # customized bar plot
#' plot(fit1, type = "bar", top = 15, horiz = TRUE, 
#'        bw = TRUE, title = "Top 15 CSMFs", angle = 70, 
#'        err_width = .5, err_size = .6)
#' 
#' # run with sub-populations
#' fit2<- insilico( data, subpop = subpop, 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE, 
#'				 auto.length = FALSE)
#' # default plot
#' plot(fit2, type = "compare", top = 5, title = "Top 5 causes comparison")
#' # customized single sub-population plot
#' plot(fit2, which.sub = "Unknown", 
#' 	 title = "Top 10 causes in HIV status unknown population")
#' # customized plot with specified causes, with abbreviation here.
#' some_causes <- c("HIV", "Pulmonary", 
#' 				"Other and unspecified infect dis")
#' plot(fit2, type = "compare", horiz = FALSE,  causelist = some_causes,
#' 		   title = "Infectious diseases in three sub-populations", 
#' 		   angle = 20)
#' }
#' 
#' @export plot.insilico
plot.insilico <- function(x, type = c("errorbar", "bar", "compare")[1], 
	top = 10, causelist = NULL, which.sub = NULL, 
	xlab = "Causes", ylab = "CSMF", title = "Top CSMF Distribution", 
	horiz = TRUE, angle = 60, fill = "lightblue", 
	err_width = .4, err_size = .6, point_size = 2, 
	border = "black", bw = FALSE, ...){
	
	sx <- summary(x)
	# sx2 <- summary(x,  CI.csmf = 0.5)

	# to please R CMD Check with ggplot
	Group <- Mean <- Lower <- Upper <- Causes <- NULL

	if(type == "compare" && is.null(sx$subpop_counts)){
		type <- "bar"
		warning("No group detected, switch to simple bar plot")
	}
	##
	## for multiple populations plot 
	##
	if(type == "compare"){
		# get which causes to plot
		# if causelist is provided then fine, 
		# otherwise construct from top argument
		if(top > 0 && is.null(causelist)){
			for(i in 1:length(sx$csmf.ordered)){
				causelist <- c(causelist, 
						setdiff(rownames(sx$csmf.ordered[[i]])[1:round(top)], 
								causelist))
			}
		}
		
		csmf.toplot <- NULL
		csmf.sub <- NULL
		for(i in 1:length(sx$csmf.ordered)){
			vcauses <- rownames(sx$csmf.ordered[[i]])
			if(class(causelist) == "character"){
				toplot <- rep(1, length(causelist))
				for(j in 1:length(causelist)){
					full <- match.arg(tolower(causelist[j]), tolower(vcauses))
					toplot[j] <- match(full, tolower(vcauses))
				}
			}else{
				# if numeric cause list is given
				toplot <- match(rownames(sx$csmf[[1]])[causelist], 
								vcauses)
			}
			csmf.toplot.tmp <- sx$csmf.ordered[[i]][toplot, c(1, 3, 5)]
			csmf.toplot <- rbind(csmf.toplot, csmf.toplot.tmp)
			csmf.sub <- c(csmf.sub,  rep(names(sx$csmf.ordered)[i], 
									     dim(csmf.toplot.tmp)[1]))
		}
		rownames(csmf.toplot) <- NULL
		csmf.toplot <- data.frame(Causes=causelist, csmf.toplot, Group = factor(csmf.sub, levels = unique(csmf.sub)))
	}
	##
	## for single population plot 
	##
	if(type == "bar" || type == "errorbar"){
		# check which subpopulation to plot first.
		# if no subpopulation
		if(is.null(sx$subpop_counts)){
			csmf.ordered <- sx$csmf.ordered[, c(1, 3, 5)]
		}else{
			# if there are multiple subpopulation
			if(is.null(which.sub)) stop("More than one groups detected. Please specify which to plot")
			if(class(which.sub) == "character"){
				which.sub <- match(which.sub, names(sx$csmf.ordered))
				if(is.na(which.sub)) stop("Invalid 'which.sub', no match found.")
			}
			csmf.ordered <- sx$csmf.ordered[[which.sub]][, c(1, 3, 5)]	
		}
		# get cause names
		vcauses <- rownames(csmf.ordered)
		# get which causes to plot
		if(is.null(causelist)){
			if(top > 0){
				toplot <- 1:round(top)
			}else{
				stop("Invalid causes to plot. Please specify causelist or a positive top")
			}
		}else{
			if(class(causelist) == "character"){
				toplot <- rep(1, length(causelist))
				for(i in 1:length(causelist)){
					full <- match.arg(tolower(causelist[i]), tolower(vcauses))
					toplot[i] <- match(full, tolower(vcauses))
				}
			}else{
				# if numeric cause list is given
				toplot <- match(rownames(sx$csmf)[causelist], 
								vcauses)
			}
		}
		csmf.toplot <- data.frame(Causes = rownames(csmf.ordered[toplot, ]), 
								  csmf.ordered[toplot, ])
	}

	
	# making bar plot
	# require: number of top causes or which causes to plot;
	#		   color list
	if(type == "bar"){
		# initialize ggplot, force order of bars
		if(horiz){
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
									 y=Mean))	
		}else{
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
									 y=Mean))
		}
		g <- g + geom_bar( stat="identity", 
				 colour=border, fill=fill,  size = .3)
		## todo		
		g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper), size = err_size, width = err_width,  position = position_dodge(.9))
		g <- g + xlab(xlab) + ylab(ylab) 
		g <- g + ggtitle(title)
		if(horiz) g <- g + coord_flip()
		if(bw) g <- g + theme_bw()
		if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))	
		return(g)
	}

	# making error bar plot
	# require: number of top causes or which causes to plot;
	#		   color list
	if(type == "errorbar"){
		# initialize ggplot, force order of bars
		if(horiz){
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
									 y=Mean))	
		}else{
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
									 y=Mean))
		}
		g <- g + geom_point( stat="identity", 
				 colour=border, fill=fill,  size = point_size)
		g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper), size = err_size, width = err_width, position = position_dodge(.9))
		g <- g + xlab(xlab) + ylab(ylab) 
		g <- g + ggtitle(title)
		if(horiz) g <- g + coord_flip()
		if(bw) g <- g + theme_bw()
		if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))	
		return(g)
	}

	# making bar plot for comparison
	if(type == "compare"){
		# initialize ggplot, force order of bars
		if(horiz){
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
									 y=Mean, ymax=max(Upper)*1.05,
									 fill = Group))
		}else{
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
									 y=Mean, ymax=max(Upper)*1.05, 
									 fill = Group))
		}
		g <- g + geom_point( aes(color=Group), position=position_dodge(0.5), size = point_size) 
		g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Group), size = err_size, width = err_width, position = position_dodge(.5))
		g <- g + xlab(xlab) + ylab(ylab) 
		g <- g + ggtitle(title)
		g <- g + scale_y_continuous() 
		if(horiz) g <- g + coord_flip()
		if(bw) g <- g + theme_bw()
		if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))
		return(g)
	}   

}
