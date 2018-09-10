%
% PREDICTIVE FLOW ANALYSIS
%
% PickLambda.m
% Purpose: Pick the optimum lambda after running CrossValidation
%
% INPUT: Summary of CVlasso fit info for lambda paramaters 
%
% OUTPUT: Save optimum lambda parameter to bestlambda.mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lasso fit info has all the necessary parameters
load('../analysis_files/Parametertesting_80subjects.mat')
% loads the CVlassoFitinfo structure 

% Read the lambda values
lambdas=CVlassoFitinfo(1,1).lambdas;

for i=1:626 % loop over ROIS
    for l=1:25 % 25 lambda values
        mse=0;
        Sum_mse=0;
        for s=1:80 %80 subjects
          
            if~isempty(CVlassoFitinfo(s,i).cvm)
            mse=CVlassoFitinfo(s,i).cvm(l,1);
            rtest=CVlassoFitinfo(s,i).Rtest(l,1);
            rtrain=CVlassoFitinfo(s,i).Rtrain(l,1);
            Sum_mse=Sum_mse+mse;
            
            rtests(l,s,i)=rtest;
            rtrains(l,s,i)=rtrain;
            else
                mse= NaN;
                Sum_mse=Sum_mse;
                rtest=NaN;
                rtests(l,s,i)=0;
                rtrains(l,s,i)=0;
            end
        end
        MSE(l,i)=Sum_mse;
    end
    [M,I] = min(MSE(:,i)); % pick the one with min Rsquare training 
    ROI_I(i,1)= I; % to check best lambda for individual ROIs
    ROI_lambda(i,1)= lambdas(I);
end

overallMSE=mean(MSE');
[bestlambdaMSE,Ibest]=min(overallMSE); % Best overall lambda 

bestlambda=lambdas(Ibest);

save('../analysis_files/bestlambda.mat','bestlambda','-v7.3')