function [D] = windowed_freq_distances(X_f,include_zero_freq)
%Compute pairwise distances between signals
% X_f: freq representation of average signal over each window


[N,~,num_windows] = size(X_f);

D = zeros(num_windows, num_windows);


D = inputArg1;
outputArg2 = inputArg2;
end

