# HBT_Analysis

This is my HBT framework to calculate the HBT correlation function. 

I have used a hybrid framework, involving MUSIC + iSS + UrQMD, to generate the particle data set. For some other model or output dataset, the files in /src/ can be modified accordingly. 

For the time being, this HBT_Analysis only computes the 1-D and 3-D correlation functions, and provides output datafiles.
* For 1-D correlation function, the output files have two columns: Qinv and C(Qinv)
* For 3-D correlation function, the output files have four columns: Qout, Qside, Qlong and C(Q)

These output files can be used to plot the correlation functions or to fit them to obtain the HBT radii values. 

I will further modify the existing codes as well as add new codes to incorporate the fitting of the correlation functions and calculate the HBT radii values. 

The repository will soon be updated. Stay tuned! 

 
