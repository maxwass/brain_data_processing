function [covs, corrs] = windowed_fcs(dtseries, windowsize, movesize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dtseries_size = size(dtseries);
num_roi   = dtseries_size(1);
num_obsvs = dtseries_size(2);

%how many covariances will we compute?
wind_end = windowsize;
num_covs = 0;
while wind_end < num_obsvs
    num_covs = num_covs + 1;
    wind_end = wind_end + movesize;
end
fprintf('   computing %d window covs of size %d, shifting by %d\n', num_covs, windowsize, movesize);

%compute covs
covs  = zeros(num_roi, num_roi, num_covs);
corrs = zeros(num_roi, num_roi, num_covs);

window_start_idx = 1;
for i = 1:num_covs
    %window is windowsize elements in covariance computation
    window_end_idx = window_start_idx + windowsize;
    cov_           = cov(dtseries(:,window_start_idx:window_end_idx)');
    covs(:,:,i)    = cov_;
    corrs(:,:,i)   = corrcov(cov_);
    
    %shift start of window over by movesize
    window_start_idx = window_start_idx + movesize;
end

%{
%calculate time series covariance trajectory
wind_start_idx = 1;
wind_end_idx   = wind_start_idx + windowsize;

idx = 0; 
while(wind_end_idx < nTR_LR)
    idx = idx + 1;
    cov_traj_LR(idx,:,:) = cov(dtseries(:,wind_start_idx:wind_end_idx)');
    wind_start_idx = wind_start_idx + movesize;
    wind_end_idx = wind_end_idx + movesize;
end
%}

end

