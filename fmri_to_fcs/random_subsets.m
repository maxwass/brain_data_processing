function [subsets] = random_subsets(signals, num_subsets, sps)
% returns 3D tensor where each slice are s randomly sampled signals
% sps :: int - number of signals to sample for each subset (samples per
% subset)
%
% num_subsets :: int - number of times to sample

% exclude thhis for now. All samples taken WITH replacement (bc more samples than observations)
%   replacement :: logical - whether to sample with replacement or not

if length(size(signals))~=2
    error('signals must be a 2D matrix, where columns are vector observations');
end

[num_roi, num_obsvs] = size(signals);

replace = true;
idxs = randsample(num_obsvs, num_subsets*sps, replace);

subsets = signals(:, idxs);

%3D tensor of stacked matrices. Each matrix slice has s (column) vector observations
%in it.
subsets = reshape(subsets, num_roi, sps, num_subsets);

%{
%matrix i <->  signals in region window i
subsets = zeros(num_roi, s, num_subsets);

for i = 1:num_subsets
    if replacement %sample with replacement
        idxs = randsample(num_obsvs, s, replacement);
    else           %sample w/o replacement
        idxs = randsample(num_obsvs, s);
    end

    subset         = signals(:,idxs);
    subsets(:,:,i) = subset;
end
%}
