Version 0.51

Removed workspace files from repo.

Verison 0.5

The code to produce the results for each drug can be found in their respective folders. All the code make use of common functions, which can be found in the "Common" folder. To reproduce the main figures, one can run the code in the "Main Figure Scripts" folder. 

For each drug, the script folder contains the code that correspond to a step in the workflow to reproduce the results.  

Bortezomib:
1) bortezomibProcess.R: preprocesses the clinical datasets, also applies ComBat

All drugs:
2) [drug name]_preparing_data.R: preprocesses the CGP dataset, homogenizes CGP and clincial dataset, and generates the random partitions 
3) pp.R, cp.R, cpp.R: runs the experiments for P2P, C2P, and CP2P
4) cp_plots.R, cpp_plots.R, pp_plots.R, ... _var.R: produces the plots for each case

Workflow:


FAQ:
Why do errors occur? 
The path is incorrect. This verison of the code contains some absolute paths that worked fine on my computer but will not work on yours. You'll have to change it to yours. 

Why is it taking so long to run? 
Results were produced on clusters of computers. On a local machine you will only realistically be able to do a few runs. Be aware that mRMR takes a lot of itme to run. 

