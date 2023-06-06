# 2022-ConstrainedLoadControlStudy
Purpose:

This repository contains MATLAB code and bash scripts used to perform the design studies in [1].
The data produced have also been added (not added yet).
These purpose of these design studies is to characterize the potential performance of a wave energy converter (WEC) given a set of constraints on the load by the power take-off (PTO).
Model predictive control is used to provide a constrainted nonlinear optimal load control.

Design studies:

The three design studies reported in [1] are exectued by running calling the following function m-files found in the root directory:
- study_MPLSparameters.m
- study_loadScheduleConstraints.m
- study_yearlyAve.m

These functions take an arguement for a specific set of parameters within the respective design study.
BASH scripts are included in teh root directory with the same name (except with extension .sh instead of .m) which execute a single iteration as a SLURM job with the intention of being used as part of a job array. 
The command to submit the job array is included as a comment with the BASH scripts themselves.
Data files produced by the scripts are saved to the root directory.
It is recomended that the root directory be clean before running the scripts and to only run one study per clean project directory.

Models/Utilities:

Functions specific to implimenting model predictive control are included in the root directory.
Functions and data files implimenting the system model are included in the 'Modeling' directory. 

Data:

Data produce by the design studies are saved as mat-files in the 'Data' directory (not added yet).
Each contains the results of the simulation in the structure 'out' and the parameters used by each simulation in the structure 'par'.

Post-processing:

Scripts used to post-process the data mat-files and produce the figures in [1] are included in the root directory. 
These should be run with the directories containing the scripts and the data mat-files added to the path. 
  
References:

[1] Simmons and Van de Ven, "Limits on the Range and Rate of Change in Power Take-Off Load in Ocean Wave Energy Conversion: A Study Using Model Predictive Control". <JournalName> Under Review <correct paper citation after acceptance>
