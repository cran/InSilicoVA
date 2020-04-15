.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'InSilicoVA' package as:\n",
         "Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C. Crampin, Kathleen Kahn and Samuel J. Clark (2016). ",
         "Probabilistic cause-of-death assignment using verbal autopsies, ",
         "Journal of the American Statistical Association, 111(515):1036-1049, DOI: 10.1080/01621459.2016.1152191.\n"),
      domain = NULL,  appendLF = TRUE )
}