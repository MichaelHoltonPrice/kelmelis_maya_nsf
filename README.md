# Overview
This README describers the steps to reproduce the results for the following
NSF proposal by Kelmelis et al.:

A new end-to-end Bayesian and osteological approach to reconstructing
demographic changes for the ancient Maya

This README assumes use of a terminal/command window, but most of the steps
remain the same if a different approach is used.

For a detailed description of the statistical methodology, see the top of
the file do_power_anaylsis.R. Additional comments and documentation are also
provided in the remainder of the code.

# Setup
Clone this repository, enter the newly created directory, and start R:

```console
https://github.com/MichaelHoltonPrice/kelmelis_maya_nsf
cd kelmelis_maya_nsf
R
```

Install required packages:

```console
install.packages('devtools')
install.packages(c('matrixStats', 'HDInterval', 'doParallel'))
devtools::install_githb('eehh-stanford/baydem')
```

# Create the simulated results
Run the following script to create simulated results.

```console
source('do_power_analysis.R')
```

This will write an output .rds file to the /data directory for each
experiment/simulation in the power calculation. It will also create a PDF plot
at the root of the repository that beings power_curve_ and includes the analysis
name in the filename. The analysis name is automatically generated from the
date and time when the script was started.