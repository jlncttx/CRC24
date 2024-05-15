Codes associated with Cattiaux, Ribes and Cariou (2024), How extreme were daily global temperatures in 2023 and early 2024?, submitted to NCC. Associated data is available here: https://zenodo.org/records/11198731

Details:

1. Running the scanning procedure.

scan_tglo.R is the executable script that gets the ERA5 data, performs the analysis and saves the results. Results are stored as "scan_1d" objects, i.e. a list containing all information of a scan with one method (data, trend, fits, p1, etc.).

scan_source.R is sourced by scan_tglo.R. It contains generic information for the scan: required packages, input / output directories, and a few useful function (e.g. for the computation of p1).

rbase.R is sourced by scan_source.R. It contains basic functions, e.g. the treatment of netcdf files (function myno() imports a netcdf in a list, etc.) or the treatment of time/dates, etc. It also loads several packages that the user will need to install.

forced.response.R is sourced by scan_source.R. It contains functions used for the detrending procedure. See Rigal et al. (2019) and Ribes et al. (2022) for further details.

scan_1d_exe.R is a subroutine of scan_tglo.R. It performs the scan at one location for one method and returns a scan_1d object.

scan_1d_sub.R is sourced by scan_tglo.R. It contains functions useful to scan one time series, eg. compute.stat.1d() that generates the nday x nduration matrix of p1, p0, etc. values.

scan_methods.txt is read by scan_1d_exe.R. It contains the namelist of parameters for each scanning method.

2. Visualizing the results.

scan_functions.R contains functions to get scan.1d objects and select the most extreme events (full chronology).

scan_figures.R contains functions to plot the results.

CRC24.R is an executable script that reproduces tables and figures of the paper.

CRC24_figures.R contains functions to generate additional figures for the paper.
