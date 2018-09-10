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


load('AAL626_final_HCP842_1mm.fib.gz.prob_connectivity.mat')
% loads the struct 'out'
% out.conn_mat field is a 626x626x500 matrix that keeps the number of streamlines
% per each ROI pair for 500 bootstrapped tractography parameter sets

% Run one-sided t-test on the distribution of streamlines for each ROI pair

% initialize matrices 
hmatrix=zeros(626,626);
tmatrix=zeros(626,626);
pmatrix=zeros(626,626);

for i=1:626 % loop thru the bootstrapped tractography matrix
    for ii=i+1:626
        [h,p,ci,stats] = ttest(squeeze(out.conn_mat(i,ii,:)));
        
        hmatrix(i,ii)=h;  %Hypothesis testing hmatrix=1 if null hypo is zero
        tmatrix(i,ii)=stats.tstat;
        pmatrix(i,ii)=p; % p value of each t-test 
        
    end
end

% Calculate the new significance level, adjust for multiple comparisons
% Using FDR.m, what is the false discovery rate for p=0.025

% First extract the upper triangle of pmatrix (obtain unique #comparisons)
mask = triu(true(size(pmatrix)),+1); % mask for upper triangle 
pvec = pmatrix(mask);  % upper triangle p vector

%Apply the fdr.m
[pID,pN]=fdr(pvec,0.025);

%fdr.m calculates two new significance threshold for different assumptions
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold
% We pick pID as the FDR thresholded significance level


% Create binary connectivity matrix
% Set to 1, if the tractography ttest is less than the adjusted significance pID

BinaryConn=zeros(626,626);
for i=1:626
    for ii=i+1:626
        if(~isnan(pmatrix(i,ii)) && lt(pmatrix(i,ii),pID))
            % check pmatrix is element not NAN and 
            % check if the pmatrix value is less than the significance threshold
            BinaryConn(i,ii)=1;
            BinaryConn(ii,i)=1;
        end
    end
end


% Sorting the Binary Connectivity matrix into left and right hemispheres
load('sort_index.mat') % 626x1 vector for sort index 
sortedBinaryConn=sortmatrix(BinaryConn,sort_index);

save('../../analysis_files/StructuralMatrix.mat','BinaryConn','-v7.3')
save('../../analysis_files/SortedStructuralMatrix.mat','sortedBinaryConn','-v7.3')




