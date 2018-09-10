# README.md 
# Analysis Pipeline for Predictive Flow Functional Connectivity Analysis

cite Fatma Uyar et. al. 2018

=========================================================================

1. Predictive flow code has the directory structure below
    
   /analysis_files  % to hold the input/output .mat files
   /glmnet_matlab % copy of cvglmnet toolbox 
   /model         % contains the Predictive Flow connectivity matrix and atlas etc
   /fmri_data      %HCP900 resting state data
   		-/preprocessed_holdout   %holdout set
   		-/preprocessed           %training and testing set
   /scripts       % matlab scripts for Predictive flow analysis 
		-/FC  % scripts for building the Predictive Flow Functional Connectivity
   		-/SC  % code for building the Structural Connectivity matrix

2. Add cvglmnet matlab toolbox to matlab path for Lasso Calculations
   'addpath(genpath('glmnet_matlab'))'
   cvglmnet toolbox reference link
        https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html

3. Run the /FC Predictive Flow Analysis scripts in the order below: 

	1) LambdaParameterCrossValidation.m 
	2) Picklambda.m
	3) PredictiveFlow_IndividualSubs.m



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% LambdaParameterCrossValidation.m
% Purpose: Run cross validated cvglmnet lasso to search for best lambda
%
% INPUT: Preprocessed holdout set subjects 
%        ROI bold signals TASK1 run1 and run2 only
%
% OUTPUT: Summary of CVlasso fit info for lambda paramaters per ROI/subject
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PickLambda.m
% Purpose: Pick the optimum lambda after running CrossValidation
%
% INPUT: Summary of CVlasso fit info for lambda paramaters 
%
% OUTPUT: Save optimum lambda parameter to bestlambda.mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PredictiveFlow_IndividualSubs.m
% Purpose: Calculate lasso beta parameters for training/test subjects 
%
% INPUT: Preprocessed training/test ROI data
%        Optimum lambda value from cross validation analysis 
%
% OUTPUT: LassoFit info per subjects per roi
%         Predictive Flow Beta Matrix( averaged over subjects)
%
% Calls for sortmatrix.m to Organize BetaMatrix into hemispheres
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



4. scripts/SC holds the code for building the structural matrix 

%
% STRUCTURAL CONNECTIVITY ANALYSIS
%
% CreateStructuralConnectivityMatrix.m script
% Purpose: Run t-test on the permuted tractography results on to create 
%          binary structural matrix 
%
% INPUT: Tractography calculations from dMRI analysis
%       "AAL626_final_HCP842_1mm.fib.gz.prob_connectivity.mat"
%        "sort_index.mat" sort ROI into hemispheric ordering
%
% OUTPUT: Binary structural connectivity matrix
% 
% Calls for fdr.m and sortmatrix.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

