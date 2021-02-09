function [subsets] = random_subsets(dtseries, num_subsets, s, replacement)
% returns 3D tensor where each slice are s randomly sampled signals
% s :: int - number of signals to sample for each subset
% num_subsets :: int - number of times to sample
% replacement :: logical - whether to sample with replacement or not

dtseries_size = size(dtseries);
num_roi   = dtseries_size(1);
num_obsvs = dtseries_size(2);

%matrix i <->  signals in region window i
subsets = zeros(num_roi, s, num_subsets);

for i = 1:num_subsets
    if replacement %sample with replacement
        idxs = randsample(num_obsvs, s, replacement);
    else           %sample w/o replacement
        idxs = randsample(num_obsvs, s);
    end

    subset         = dtseries(:,idxs);
    subsets(:,:,i) = subset;
end

