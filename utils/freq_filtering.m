function [filter_signals_freq] = freq_filtering(signals_freq,freq_intervals)
%Given signals in frequency domain (columns are signals), take 
%components from non-intersecting intervals in frequency domain to create
%fitlered signal

%signals_freq :: Nxnum_obsvs [float]
%intervals :: kx1 *cell array* of 1x2 [int]'s


%% error checking
if length(size(signals_freq))~=2
    error('requires 2D array where columns are observations');
end
[N, ~] = size(signals_freq);

last_high = -1;
for j = 1:length(freq_intervals)
    interval = freq_intervals{j};
    low = interval(1);
    high = interval(2);
    
    if low>high %has to be valid interval
        error('invalid %dth interval given in freq_intervals: %d > %d', ...
            j, low, high);
    end
    if low <= last_high %should be non-oversecting with previous interval
        error('intersecting interval given between intervals %d & %i: %d <= %d', ...
            j-1, j, low, last_high);
    end
end
    
%ensure highest interval is not out of bounds
last_interval = freq_intervals{end};
if last_interval(2)>N
    error('freq range too high (%d): only %d eigenvalues', last_interval(2), length(eigs));
end


%% construct filtered signal
filter_signals_freq = zeros(size(signals_freq));
            
for i = 1:length(freq_intervals)         
	interval = freq_intervals{i};
    [low,high] = deal(interval(1), interval(2));
    
    %place freq components between low and high into fitlered signals
    filter_signals_freq(low:high,:) = signals_freq(low:high,:);     
end


end

