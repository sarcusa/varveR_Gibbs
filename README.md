Supporting code for New approaches to dating intermittently varved sediment, Columbine lake, Colorado, USA by S. Arcusa, N.P. McKay, C. Wiman, S. Patterson, S.E. Munoz and M.A. Aquino-Lopez.

The following files contain functions which should be saved in a separate folder:
GibbsrRelatedFunctions.R, ImportFiles.R, plotting.R, simulateOverAndUnderCounting.R, varveModel.R

Code is to be run in the following order

1. script_1 contains code to run the first half of the analysis with asymmetrical and symmetrical priors. It prepares the data for the Gibbs sampler.
2. Gibbs_slurm.R contains code that uses output from script_1 to produce the datafiles for script_2. Will need to change the directory paths. May need to create a launch file and adjust for specific super computer 
  Best to run Gibbs_slurm.R at least 50,000 times.
3. script_2 will combine observers and cores.

Other code files, to be run in this order, for extra analyzes and plots.

1. COL_script.R contains code to create the age-depth model plots.
2. multi_cores_observers.R uses the output of COL_script.R to prepare the data frames of laminations from multiple cores and observers. 
3. Use Gibbs_sampling_model.R to start the gibbs model. Run once for each observer. Only run for a few hundred or a thousand iterations max. Use Gibbs_slurm.R for thousands of iterations.





