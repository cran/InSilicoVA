# InSilicoVA - changes

Version 1.2.7 (2020-02-18)
==========================
+ Include system checks for Java version at least 7.0.

Version 1.2.6 (2020-02-17, on CRAN)
==========================
+ Include function to load the latest P(S|C) from InterVA-5 software

Version 1.2.6 (2018-12-10, on GitHub)
==========================
+ Bug fix for model with customized P(S|C)
+ Bug fix for scenario where all deaths due to external causes in WHO 2016 format.

Version 1.2.5 (2018-10-28)
==========================
+ Bug fix for subsetting data with (1) no external causes, and (2) only external causes.
+ Added checks for extra columns in physician debias function.

Version 1.2.4 (2018-09-12)
==========================
+ Add option to exclude impossible causes based on different rules.
+ The previous default rule to exclude impossible causes is now called 'subset'. It further include impossible causes screening based on the following symptoms: 
    + baby died within 24 hours; between 24 and 48 hours; between 48 hours and within one week; between one to four weeks. 
    + baby born dead.
    + woman 12-19 years old; 20-34 years old; 35-49 years old.
+ Bug fix symptoms not negated when absence is taken as the significant value. 
+ Bug fix that converts absence to missing for WHO 2016 input introduced in version 1.2.3.
+ Bug fix for individual probability calculation when no cause has positive probabilities.


Version 1.2.3 (2018-08-27)
==========================
+ Add codes to safely handle data containing NA values.
+ Fix the writing of warning and error log files for both WHO 2012 and 2016 formats, and allowed customized directory. 
 

Version 1.2.2 (2018-07-13)
========================== 
+ Fix typo in checking WHO 2016 input using InterVA5 rules, and update the check rules for external causes based on all symptoms associated with external causes (instead of using only the most deterministic one as in the previous versions).
+ Fix typo in shifting columns for WHO 2016 input.
+ Fix the issue that sometimes 0 acceptance rate is displayed when the chain is in fact moving.
+ Fix inconsistent names in InterVA5 default input.
+ Removed most of the cat() calls and replaced them with message().
+ Add subpop selection option for 'compare' plots.


Version 1.2.1 (2018-07-06)
==========================
+ Fix typo in processing WHO 2016 input.
+ Fix inconsistency from symptom name change (sk\_les and skin\_les).

Version 1.2.0 (2018-04-20)
========================== 
+ Add WHO 2016 support using InterVA5 input
+ Add option to allow group code in output
+ Add warning and error logs in the InSilico object, so that removed observations can be traced back
+ Make probbase consistent with InterVA-4.03 instead of 4.02. Changes should be minimum, but good to be consistent.
+ Allow user defined impossible causes to be removed


Version 1.1.5 (2017-11-17)
==========================
+ Fix issue for data without external causes.
+ Fix issue with print methods after last update.
+ Fix issue with neonate and child death assigned to impossible external causes.
+ Change data check steps to update symptoms to be missing instead of absence to introduce symptom dependence structures. 
+ Automatic remove impossible causes from CSMF based on subpopulation gender and age. 


Version 1.1.4 (2017-01-24)
========================== 
+ Fix issue with sub-population specification and output labels of non-standard InterVA-4 input.
 

Version 1.1.3 (2017-01-02)
==========================
+ Fix issue with label order in probbase output.
+ Allow input "Yes" to be denoted with either "y" or "Y".
+ Default data checking rules to be consistent with InterVA-4.03 instead of the previous version.
+ Fix typo in document.

    
