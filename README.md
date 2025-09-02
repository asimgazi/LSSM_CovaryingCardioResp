# LSSM_CovaryingCardioResp

This repository contains code to replicate the methods described in the below referenced paper.

A. H. Gazi et al., "Modeling Latent Dynamics of the Autonomic Nervous System in Response to Trauma Recall and Non-Invasive Vagus Nerve Stimulation," in *IEEE Transactions on Biomedical Engineering*, doi: [10.1109/TBME.2025.3580051](https://doi.org/10.1109/TBME.2025.3580051).

**The primary purpose of this repository is to provide reference code to replicate methods for other applications.** 

If you would also like to see how this specific code runs, the `main.mlx` script is the place to start. The `main.mlx` script is a MATLAB live script - similar to a Python notebook. Below are general steps to run the code, noting that aspects of the code will need to be modified based upon the data you have available.

1. Make sure you have all the necessary MATLAB toolboxes: system identification toolbox, signal processing toolbox, statistics and machine learning toolbox, and the parallel computing toolbox if you would like to stick with `parfor` loops (if not, simply change `parfor` loops to `for` loops and replace `parsave` with your desired saving function)
2. Update directories, variable names, etc. in the script based upon where you have processed physiological features and protocol timings saved, where you would like results to be saved, how you would like results to be named, etc.
3. Run the script and store any necessary results for further statistical analyses
4. Run the additional statistical analyses scripts (in R)

Note that much more analysis was conducted and many more plots were generated compared to the figures and analyses included in the paper. If you only want to generate the results shown in the paper, identify those specific code blocks and run those. 
