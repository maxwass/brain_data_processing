% create fc window datasets

% use parallel for loop

%% Pipeline
% load data

%create fc's w/ optional filter
% create average signals in window
[signal_windows] = windowed_signals(dtseries, windowsize, movesize);

% NEED
% filter
%[metrics] = some metric of ave_signals

% pre ave filters


% post ave filters
% ex) only take ave signals with energy less than average
%        energies = norm of cols 
%        above_ave_energies = energies >= mean(energies)
%        energies(above_ave_energies) = []
% ex) low pass filter ave signals
%        freq ave signals = GFT(ave_signals)
%        lpf_freq_ave_signals = set high freqs to zero
    

% compute covariances


%save to individual mat file