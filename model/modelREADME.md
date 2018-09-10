# modelREADME.md 
# data files from the Predictive Flow Functional Connectivity Analysis

cite Fatma Uyar et. al. 2018

=========================================================================

LIST OF FILES IN THE model FOLDER 

Binary structural connectivity matrix for AAL626 regions: SortedStructuralMatrix.mat
 
Predictive Flow Functional connectivity matrix ( Average LASSO parameter matrix) 
SortedLassoBetaMatrix.mat

The vector of mean r-square values for all 626 regions: meanlasRsq.mat

The vector of standard error r-square values for all 626 regions: stderrlasRsq.mat

The spreadsheet with the AAL labels, their MNI locations, and sorted ID ordering: AAL626_atlas.xlsx


NOTES ON PREDICTIVE FLOW and model files
the LASSO parameter matrix is a measure of functional connectivity that doesn't use correlation. 

Each row of the matrix is an estimated regression model for a specific ROI, where the off diagonal entries are the fitted regression coefficients from a LASSO regression. 

So for each subject we run 626 independent regression models where we learn the association between all n-1 nodes (independent variables) and the timeseries at the target node (dependent variable). 

The lasso matrix is the average across 700+ HCP subjects. 

We learned the LASSO model from one resting state fMRI run and evaluated its accuracy on the second resting state fMRI run with r-squared measure. 

Thus, the r-squared vectors show the predicted vs. observed accuracies of each ROI from a hold out test set. 