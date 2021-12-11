rAmpsMsre
=========

[![GitHub release](https://img.shields.io/github/release/Alexeika2/rAmpsMsre.svg)](https://github.com/Alexeika2/rAmpsMsre/releases)

# Introduction

*rAmpsMsre* is an [R](http://en.wikipedia.org/wiki/R_%28programming_language%29) project
for DNA methylation analysis and select virtual amplicons for MSRE-PCR based on
hierachical clustering with custom distance function between MSRE recognition sites
including manhattan distance on methylation in sample and physical distance in genome.
Plots designed to assess amplicon primer design with visual control.

The package is designed to deal with
[bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) output files and
[rrbsData](https://github.com/tanas80/rrbsData) S4 structures.

## Current Features

 * Loads methylation calling data based on data.frame sample description
 * Easy R-style filtering object with indexes of specimens or CpG-dinucleotides
 * Two layered matrix data design (B-value and coverage).
 * Some reports

## Installation
In R console:
```r
install.packages("BiocManager")
library(devtools)
install_github("tanas80/rrbsData", build_vignettes=FALSE,
               dependencies=TRUE, ref="init", repos=BiocManager::repositories())
install_github("Alexeika2/rAmpsMsre", build_vignettes=FALSE,
               dependencies=TRUE, ref="init", repos=BiocManager::repositories()))

```
-------
# How to Use

[example_amps.R](https://github.com/Alexeika2/rAmpsMsre/blob/master/examples/example_amps.R)
