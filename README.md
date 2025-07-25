# HBT_Analysis

This is my HBT framework to calculate the HBT correlation function. 

I have used a hybrid framework, involving MUSIC + iSS + UrQMD, to generate the particle data set. For some other model or output dataset, the files in `src/` can be modified accordingly. 

For the time being, this HBT_Analysis only computes the 1-D and 3-D correlation functions, and provides output datafiles.
* For 1-D correlation function, the output files have two columns: Qinv and C(Qinv)
* For 3-D correlation function, the output files have four columns: Qout, Qside, Qlong and C(Q)

These output files can be used to plot the correlation functions or to fit them to obtain the HBT radii values. 

## Repository Overview

Below, I provide an overview of the different folders present in this repository, as well as what each folder contains:

| Folder           | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| `config/`        | Contains config.ini file for providing parameter configuration values.|
|`input_data/`     | Include the particle data set from which all necessary data are to be extracted.|
| `sample_outputs/`| Contains generated results, including graphs, tables, and processed data. These are mere samples used to verify proper working of codes and were not used for original analysis purposes.|
| `scripts/`       | Contains bash script as a separate way to execute the entire framework without using Makefile.|
| `src/`           | Core source codes, header files and main.cpp file for analysis, computations, and generating output files.|
|`HBT_Fit\`        | Includes both Python and ROOT files to extract data from the correlation function outputs files and fit the data to obtain the HBT radii values. Also, plots the HBT radii with k_T.|


I will further modify the existing codes as well as add new codes to incorporate the fitting of the correlation functions and calculate the HBT radii values. 

The repository will soon be updated. Stay tuned! 

 
