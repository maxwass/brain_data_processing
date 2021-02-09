function [covs, ave_signals] = windowed_fcs(dtseries, windowsize, movesize)
%Scan window over all fc signal observations and compute average signals, 
% covariances, and correlations of each window

dtseries_size = size(dtseries);
num_roi   = dtseries_size(1);
num_obsvs = dtseries_size(2);

%how many covariances will we compute?
num_windows = idivide( int64(num_obsvs-windowsize), int64(movesize) ) + 1;
fprintf('   computing %d window covs of size %d, shifting by %d\n', num_windows, windowsize, movesize);

%compute covs
covs  = zeros(num_roi, num_roi, num_windows);

%column i <-> averaged signal in region i
ave_signals = zeros(num_roi, num_windows);

window_start_idx = 1;
window_end_idx   = window_start_idx + windowsize - 1; %want windowsize # elements
for i = 1:num_windows
    %window is windowsize elements in covariance computation
    %ave signal over this window
    signals = dtseries(:,window_start_idx:window_end_idx);
    ave_signals(:,i) = mean(signals,2);
    covs(:,:,i)  = cov(signals');
    
    %shift start of window over by movesize
    window_start_idx = window_start_idx + movesize;
    window_end_idx   = window_start_idx + windowsize - 1;
end


end

