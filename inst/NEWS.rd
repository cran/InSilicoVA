% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/Matrix/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{InSilicoVA}}
\encoding{UTF-8}

\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\sspace}{\ifelse{latex}{\out{~}}{ }}

\section{Changes in version 1.1.5 (2017-11-17)}{
 \itemize{
   \item Fix issue for data without external causes.
   \item Fix issue with print methods after last update.
   \item Fix issue with neonate and child death assigned to impossible external causes.
  \item Change data check steps to update symptoms to be missing instead of absence to introduce symptom dependence structures. 
    \item Automatic remove impossible causes from CSMF based on subpopulation gender and age. 
 }
}
\section{Changes in version 1.1.4 (2017-01-24)}{
 \itemize{
   \item Fix issue with sub-population specification and output labels of non-standard InterVA-4 input.
 }
}
\section{Changes in version 1.1.3 (2017-01-02)}{
  \itemize{ 
    \item Fix issue with label order in probbase output.
    \item Allow input "Yes" to be denoted with either "y" or "Y".
    \item Default data checking rules to be consistent with InterVA-4.03 instead of the previous version.
    \item Fix typo in document.

    }
}
