function [basis, varExp, proj, covMat] = func_pca(mat,varargin)
%func_pca - perform pca, return basis, variance explained, correlation matrix
%
% Syntax: [basis, varExp, corrMat] = func_pca(mat)
%
% Long description
mat = mat - repmat(mean(mat,2),1,size(mat,2));
covMat = cov(mat');
[u,v] = eig(covMat);
v = diag(v);
[varExp,idx] = sort(v,'descend');
basis = u(:,idx);
varExp = varExp ./ sum(varExp);
cumVarExp = cumsum(varExp); 
proj = basis' * mat;

end