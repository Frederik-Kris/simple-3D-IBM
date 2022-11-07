function n = get1D_index(i,j,k,NI,NJ)
%GET1D_INDEX Get the single index to access array element (i,j,k).
%   Matlab can access a multidimensional array by using multiple indices or one index.
%   This function returns the single index to use to access the same element as you get
%   if you access using the 3D indices i,j,k.
n = (k-1)*NI*NJ + (j-1)*NI + i;
end

