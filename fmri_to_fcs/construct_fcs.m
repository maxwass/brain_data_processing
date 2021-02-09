function [fcs] = construct_fcs(subsets)
%% for given tensor of subsets of (processed) fmri signals, construct fcs 
%  (covs) for each slice
%% inputs:
% subsets :: num_roi x num_samples x num_subsets double array
%   -each slice is a matrix whose columns are fmri observations
%% outputs:
% fcs = covariance matrices computed over constructed dtseries subsets

% Applies cov to each slice of tensor. Constructs new tensor of covariances.
% cov expects observations to be rows. We have columns. Transpose each
% slice before cov()
cov_t = @(x) cov(x');
fcs = apply_to_tensor_slices(cov_t, subsets);

end