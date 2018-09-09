function [sortA,i2,j] = myunique(A)
%%
% A multidimensional array
matlabversion = version;
if str2double(matlabversion(end-5:end-2)) > 2012
    [sortA, i2, j] = unique(A,'rows','legacy'); %#ok<*ASGLU>
else
    [sortA, i2, j] = unique(A,'rows');
end