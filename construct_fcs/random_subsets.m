function [subsets] = random_subsets(signals, num_subsets, sps, with_replacement)
% returns 3D tensor where each slice are sps randomly sampled signals
% sps :: int - number of signals to sample for each subset (samples per
% subset)
%
% num_subsets :: int - number of times to sample
% sps :: int - samples per slice
% with_replacement :: logical - whether or not the sampling is done with
% replacement or not


if length(size(signals))~=2
    error('signals must be a 2D matrix, where columns are vector observations');
end

[num_roi, num_obsvs] = size(signals);

%3D tensor of stacked matrices. Each matrix slice has sps (column) vector observations
%in it.
subsets = zeros(num_roi, sps, num_subsets); 
for i = 1:num_subsets
	%must resample each time b/x if sampling w/o replacement, may run out of
	%samples.
	idxs           = randsample(num_obsvs, sps, with_replacement);
	subset         = signals(:,idxs);
	subsets(:,:,i) = subset;
end
