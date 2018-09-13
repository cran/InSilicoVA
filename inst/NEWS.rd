% Check from R:
\name{NEWS}
\title{News for \R Package \pkg{InSilicoVA}}
\encoding{UTF-8}

\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\sspace}{\ifelse{latex}{\out{~}}{ }}


\section{Changes in version 1.2.4 (2018-09-12)}{
 \itemize{
   \item Add option to exclude impossible causes based on different rules.
   \item The previous default rule to exclude impossible causes is now called 'subset'. It further include impossible causes screening based on the following symptoms: 
        \itemize{
            \item baby died within 24 hours; between 24 and 48 hours; between 48 hours and within one week; between one to four weeks. 
            \item baby born dead.
            \item woman 12-19 years old; 20-34 years old; 35-49 years old.
        }
  \item Bug fix symptoms not negated when absence is taken as the significant value.      
  \item Bug fix that converts absence to missing for WHO 2016 input introduced in version 1.2.3.
  \item Bug fix for individual probability calculation when no cause has positive probabilities.
 }
}

\section{Changes in version 1.2.3 (2018-08-27)}{
 \itemize{
   \item Add codes to safely handle data containing NA values.
   \item Fix the writing of warning and error log files for both WHO 2012 and 2016 formats, and allowed customized directory. 
 }
}

\section{Changes in version 1.2.2 (2018-07-13)}{
 \itemize{
   \item Fix typo in checking WHO 2016 input using InterVA5 rules, and update the check rules for external causes based on all symptoms associated with external causes (instead of using only the most deterministic one as in the previous versions).
   \item Fix typo in shifting columns for WHO 2016 input.
   \item Fix the issue that sometimes 0 acceptance rate is displayed when the chain is in fact moving.
   \item Fix inconsistent names in InterVA5 default input.
   \item Removed most of the cat() calls and replaced them with message().
   \item Add subpop selection option for 'compare' plots.

 }
}

\section{Changes in version 1.2.1 (2018-07-06)}{
 \itemize{
   \item Fix typo in processing WHO 2016 input.
   \item Fix inconsistency from symptom name change (sk\_les and skin\_les).
 }
}

\section{Changes in version 1.2.0 (2018-04-20)}{
 \itemize{
   \item Add WHO 2016 support using InterVA5 input
   \item Add option to allow group code in output
   \item Add warning and error logs in the InSilico object, so that removed observations can be traced back
   \item Make probbase consistent with InterVA-4.03 instead of 4.02. Changes should be minimum, but good to be consistent.
   \item Allow user defined impossible causes to be removed
 }
}

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
