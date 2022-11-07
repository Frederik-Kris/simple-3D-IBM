function [i,j,k] = get3D_indices(n,NI,NJ)
%GET3D_INDICES Get the three indices to access array element n.
%   Matlab can access a multidimensional array by using multiple indices or one index.
%   This function returns the three indices to use to access the same element as you get
%   if you access using the single index 'n'.
k = floor((n-1)/(NI*NJ)) + 1;
n = n - (k-1)*NI*NJ;
j = floor((n-1)/NI) + 1;
i = n - (j-1)*NI;
end

