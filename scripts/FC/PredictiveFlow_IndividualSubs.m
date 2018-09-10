%
% PREDICTIVE FLOW ANALYSIS
%
% PredictiveFlow_IndividualSubs.m
% Purpose: Calculate lasso beta parameters for training/test subjects 
%
% INPUT: 1)Preprocessed training/test ROI data
%        2) Optimum lambda value from cross validation analysis 
%
% OUTPUT: 1) LassoFit struct that has lasso info per subjects per roi
%         2) Average beta matrix constructed from individual lasso 
%         3) StdBeta matrix calculated across subjects
%
%         LassoFit(N,626) struct has fields:
%             --betas    % individual lasso betas for each ROI
%             --Rtrain 
%             --Rtest 
%             --Y 
%             --DF 
%             --lambda 
%             --sub 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use best lambda value from cross validation set
load('../analysis_files/bestlambda.mat')
bestlam=bestlambda;

%Read in preprocessed subject ROI bold signals TASK1 run1 and run2 only
cd ../fmri_data/preprocessed_holdout
s=what;
matfiles=s.mat;
N=numel(matfiles); %number of subjects


%Define the function for rsquare calculation
rsquare = @(y, f) (1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));

% Initialize struct to keep CVglmnet calculations per lambda/ROI/subject
lassoFitinfo(N,626)=struct; 

%Start lasso calculations for predictive flow analysis
for a=1:N  %number of subjects
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
    
    opts.lambda=bestlam; % Use cross validated lambda value 
    
    for i=1:626 %626 rois
        
        traindata=subrun1;
        traindata(i,:)=0; % Set the ith ROI to 0 train on this 626*(840)
        
        testdata=subrun2;
        testdata(i,:)=0; % Remove the ith ROI to test on this 626*(840)
        
        traindata=traindata'; % each column represents one predictor 840*626
        testdata=testdata';
        
        roirun1=subrun1(i,:)'; % BOLD response of ith ROI from run1
        roirun2=subrun2(i,:)';
        
        if ~isempty(nonzeros(roirun1)) %Some of the ROI BOLDS are zero, check for that
            fit = glmnet(traindata, roirun1,[],opts);  % Run glmnet lasso
            
            ncoef=glmnetCoef(fit,bestlam);
            %this is size 627*1, first one is intercept, the others are beta
            
            predrun1 = glmnetPredict(fit,traindata,bestlam);
            predrun2 = glmnetPredict(fit,testdata,bestlam);
            
            RSq_training=rsquare(roirun1,predrun1);
            RSq_testing=rsquare(roirun2,predrun2);
            
            lassoFitinfo(a,i).betas=ncoef(2:627,:)'; % first coefficient is intercept
            lassoFitinfo(a,i).Rtrain=RSq_training;
            lassoFitinfo(a,i).Rtest=RSq_testing;
            lassoFitinfo(a,i).Y=roirun2;
            lassoFitinfo(a,i).DF=fit.df;
            lassoFitinfo(a,i).a0=fit.a0;
            lassoFitinfo(a,i).lambda=fit.lambda;
            lassoFitinfo(a,i).sub=names;

        end
 
    end
   
    
end

 save('../../analysis_files/BestLambda_IndLassoFit.mat','lassoFitinfo','-v7.3')
 
% Build a 626X626 beta matrix for each subject by aligning ROI beta rows
indBeta=zeros(626,626); % beta matrix for each individual

subjectbetamats(N,1)=struct; %struct to keep individual subject betas

for a=1:N % Looping through all subjects
    for  i=1:626 % Looping through rois
        
        indBeta(i,:)=lassoFitinfo(a,i).betas;
    end
    
    subjectbetamats(a,1).Beta=indBeta;
    
end


% Calculate the average beta matrix across all individual subject beta
% matrices
MeanBetamatrix=zeros(626,626);
stdBetamatrix=zeros(626,626);

%for each i,j element of BETA matrix,create a vector holding subject
%indeces
betavals=zeros(N,1); % create a vector 

for i=1:626
   
    for j=1:626
        for k=1:N % looping through all subjects 
         
            betavals(k,1)=subjectbetamats(k,1).Beta(i,j);
        end
        
        MeanBetamatrix(i,j)=mean(betavals);
        stdBetamatrix(i,j)=std(betavals);
    end
end



load('sort_index.mat') % 626x1 vector for sort index 
sortedBetaConn=sortmatrix(MeanBetamatrix,sort_index);

save('../../analysis_files/MeanBetaMatrix_unsorted.mat','MeanBetamatrix','-v7.3')
save('../../analysis_files/stdBetaMatrix_unsorted.mat','stdBetamatrix','-v7.3')

save('../../analysis_files/sortedBetaMatrix.mat','sortedBetaConn','-v7.3')



