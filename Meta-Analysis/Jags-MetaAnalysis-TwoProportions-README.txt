This file is Jags-MetaAnalysis-TwoProportions-README.txt

-----------------------------------------------------------

Copyright and license: Please see the accompanying file, license.txt, for current copyright and license information. As of November 17, 2016, the license.txt was as follows:

Modified MIT License

Copyright (c) 2016 John K. Kruschke

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

1. The following article is cited, with details updated appropriately: Kruschke, J. K., and Liddell, T. M. (2016, accepted pending final revision). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power analysis from a Bayesian perspective. Psychonomic Bulletin & Review. https://osf.io/dktc5/

2. The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

-----------------------------------------------------------

R scripts for Bayesian meta-analysis of success counts in control and treatment groups.

You must have the current version of JAGS installed. You must also have the current version of the runjags package installed. See 
https://sites.google.com/site/doingbayesiandataanalysis/software-installation
(Programmed in the style of the book: Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.)

The top-level script to run in R is Jags-MetaAnalysis-TwoProportions.R. In turn, that script sources Jags-MetaAnalysis-TwoProportions-Utilities.R, which R will look for in its current working directory. The utilities script will attempt to load various packages if you don't already have them.

When you run the top-level script, you can edit the code to use either of two sets of data. To load the towel-reuse data, on the line that has
dataSource = c( "TowelReuse" , "BetaBlocker" )[2]
change [2] to [1].

-----------------------------------------------------------

The R scripts were run successfully using JAGS 4.2.0, runjags 2.0.4-2, and R 3.3.0 on Windows 10: 

> sessionInfo()

R version 3.3.0 (2016-05-03)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] runjags_2.0.4-2 rjags_4-6       coda_0.18-1    

loaded via a namespace (and not attached):
[1] rsconnect_0.4.3 tools_3.3.0     grid_3.3.0      lattice_0.20-33

-----------------------------------------------------------
