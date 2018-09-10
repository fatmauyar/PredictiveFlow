%
% PREDICTIVE FLOW ANALYSIS
%
% LambdaParameterCrossValidation.m
% Purpose: Run cross validated cvglmnet lasso to search for best lambda
%
% INPUT: Preprocessed holdout set subjects
%        ROI bold signals TASK1 run1 and run2 only
%
% OUTPUT: Summary of CVlasso fit info for lambda paramaters per ROI/subject
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read in preprocessed subject ROI bold signals TASK1 run1 and run2 only
cd ../fmri_data/preprocessed_holdout;
s=what;
matfiles=s.mat;
N=numel(matfiles); %number of subjects in the holdout set


%set lambda search space
setlam=exp(linspace(log(0.0006),log(8),25));

%create folds from 840 time points
fold_id=zeros(840,1); %each fold gets 168 consecutive time points
fold_id(1:168,1)=1;
fold_id(169:336,1)=2;
fold_id(337:540,1)=3;
fold_id(541:672,1)=4;
fold_id(673:840,1)=5;

%Define the function for rsquare calculation
rsquare = @(y, f) (1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));

%Vectors to keep Rsquared of prediction for fitting to training/test data
% for each value of lambda parameter
RSq_training=zeros(size(setlam')); %n subjects 626 ROIS
RSq_testing=zeros(size(setlam'));

%initialize the matrices for training/testing data
traindata=zeros(626,840);
testdata=zeros(626,840);

% Initialize struct to keep CVglmnet calculations per lambda/ROI/subject
CVlassoFitinfo(N,626)=struct;

disp('Starting on hold out subjects cross validated lasso')
for a=1:N % for N holdout subjects
    a
    matObj=load(char(matfiles(a)));
    names= fieldnames(matObj);
    subdata = getfield(matObj,names{1});
    
    %subdata is 626*(1680),
    %first 840 columns are run1, second 840 columns are run2
    subrun1=subdata(:,1:840);
    subrun2=subdata(:,841:1680);
    
    subrun1=zscore(subrun1,0,2); % apply zscore along the rows
    subrun2=zscore(subrun2,0,2); %apply zscore along the rows
    
    
    opts.lambda=setlam;
    % Apply 5-fold lasso to get the ROI loadings
    for i=1:626 %626 rois
        
        fprintf('%d ', i);
        
        traindata=subrun1;
        traindata(i,:)=0; % zero ith ROI itself to train on this 626*(840)
        
        testdata=subrun2;
        testdata(i,:)=0; % zero ith ROI to test on this 626*(840)
        
        traindata=traindata'; % each column represents one predictor 840*626
        testdata=testdata';
        
        roirun1=subrun1(i,:)'; % BOLD response of ith ROI from run1
        roirun2=subrun2(i,:)';
        
        if ~isempty(nonzeros(roirun1)) %Some of the ROI BOLDS are zero, check for that
            cvfit_foldid = cvglmnet(traindata, roirun1,[],opts,[],5,fold_id);  %this one has fold id
            lambdas=cvfit_foldid.lambda;
            ncoef=cvglmnetCoef(cvfit_foldid,lambdas); %this is size 627*1, first one is intercept, the others are beta
            
            [lambdalen,~]=size(lambdas);
            for l=1:lambdalen
                lambda=lambdas(l,1);
                predrun1 = cvglmnetPredict(cvfit_foldid,traindata,lambda);
                predrun2 = cvglmnetPredict(cvfit_foldid,testdata,lambda);
                
                RSq_training(l,1)=rsquare(roirun1,predrun1);
                RSq_testing(l,1)=rsquare(roirun2,predrun2);
                
            end
            
            CVlassoFitinfo(a,i).betas=ncoef;
            CVlassoFitinfo(a,i).Rtrain=RSq_training;
            CVlassoFitinfo(a,i).Rtest=RSq_testing;
            CVlassoFitinfo(a,i).Y=roirun2;
            CVlassoFitinfo(a,i).DF=cvfit_foldid.glmnet_fit.df;
            CVlassoFitinfo(a,i).lambdas=cvfit_foldid.lambda;
            CVlassoFitinfo(a,i).cvm=cvfit_foldid.cvm;
            CVlassoFitinfo(a,i).sub=names;
            
        end
        
    end
   
    
end


 save('../../analysis_files/Parametertesting_80subjects.mat','CVlassoFitinfo','-v7.3')
