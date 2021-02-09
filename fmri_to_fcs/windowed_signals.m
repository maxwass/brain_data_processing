function [signal_windows] = windowed_signals(dtseries, windowsize, movesize)
% returns 3D tensor where each slice are the signals over that window


dtseries_size = size(dtseries);
num_roi   = dtseries_size(1);
num_obsvs = dtseries_size(2);

%how many covariances will we compute? 
% currently ignoring last window if it isn't full size
num_windows = idivide( int64(num_obsvs-windowsize), int64(movesize) ) + 1;
fprintf('   computing ave signals over %d windows of size %d, shifting by %d\n', num_windows, windowsize, movesize);

%matrix i <->  signals in region window i
signal_windows = zeros(num_roi, windowsize, num_windows);

window_start_idx = 1;
window_end_idx   = window_start_idx + windowsize - 1; %want windowsize # elements
for i = 1:num_windows
    %window is windowsize elements in covariance computation

    signals = dtseries(:,window_start_idx:window_end_idx);
    signal_windows(:,:,i) = signals;
    
    %shift start of window over by movesize
    window_start_idx = window_start_idx + movesize;
    window_end_idx   = window_start_idx + windowsize - 1;
end




end

