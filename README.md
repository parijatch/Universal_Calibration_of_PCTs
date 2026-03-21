# Universal calibration of heavy-tailed combination tests

This repository contains R and Julia code used to generate the simulation and application results presented in [1].

## Structure

### `simu/` contains scripts to reproduce the simulations. 

The script `run-simu.R` will generate and save simulated p-values. Then run `plot.R` to generate the plots. The script `plot-pairwise.R` generates the pairwise plots for combined p-values. 

### `nhanes/` contains code used to illustrate the Pareto combination test using NHANES data. 

This directory contains all code needed to reproduce the NHANES data analysis, illustrating the use of the Pareto Combination Test to assess the independence between multivariate health phenotypes using NHANES data.

Set the year and data directory in the file `configure.jl`.  Then run `get_data.jl` to download the files, and run `corr_proj.jl` to run the full analysis (this will take approximately one hour). Finally, run the `tables.jl` script to generate the latex output.

The raw data files and documentation are available [here](https://wwwn.cdc.gov/nchs/nhanes).

## Reference

[1] Chakraborty, P., Guo, F. R., Shedden, K., & Stoev, S. (2025). *On the universal calibration of heavy-tailed combination tests.* arXiv preprint [arXiv:2509.12066](https://arxiv.org/abs/2509.12066).
