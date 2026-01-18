# NHANES

This directory contains all code needed to reproduce the NHANES data analysis, illustrating the use
of the Pareto Combination Test to assess the independence between multivariate health phenotypes
using NHANES data.

Set the year and data directory in the file _configure.jl_.  Then run _get_data.jl_ to download the files,
and run _corr_proj.jl_ to run the full analysis (this will take approximately one hour).  Finally, run the
_tables.jl_ script to generate the latex output.

The raw data files and documentation are available [here](https://wwwn.cdc.gov/nchs/nhanes).

