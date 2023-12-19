# plmmLasso

This repository contains the implemented algorithm of a partial linear mixed-effects (PLMM) method for longitudinal data. The names and descriptions are shown below:

Code/CreateBases.R: The R script contains functions to generate the bases functions to model the nonlinear functions.

Code/CreateSim.R: The R script contains functions to simulate datasets. And a function "plot.fit" to vizualize the estimated trajectories and nonlinear functions of the estimated model with the true values.

Code/example.R: The R script demonstrates the usage of the implemented algorithm using a simulated dataset.

Code/modelselection.R: The R script contains a function that runs our method for mulitple values of the hyperparameters, compute a criteria for selection (BIC or BICC or EBIC) and return the model with the lowest selected criteria.

Code/plmmlasso.R: The R script contains the main function "plmmlasso" for model fitting.

Code/posi.R: The R script contains the function debias.plmm that allows to do post-selection inference on the fixed-effects.

Code/test.nonlinear.functions.R: The R script contains fucntions to compute bootstrapped p-values for the overall test of equality of the nonlinear functions and bootstrapped joint confidence intervals for the difference betwee the nonlinear functions.
