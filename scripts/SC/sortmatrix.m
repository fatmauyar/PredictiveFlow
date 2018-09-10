function [ sortedmatrix] = sortmatrix( matrix, sort_index )

% Resorting the ROI connectivity matrices according to sort_index of the atlas
% sort_index is (Nx1) row vector
% matrix is NxN matrix
%
Connectivity=matrix; 

[N,~]=size(sort_index); 
for i=1:N
    index=sort_index(i,1);
    sortedConnectivity(i,:)=Connectivity(index,:);%First sort rows
    
end

for i=1:N
    index=sort_index(i,1);
    ssConnectivity(:,i)=sortedConnectivity(:,index);%Then sort columns
end

sortedmatrix=ssConnectivity;
end

