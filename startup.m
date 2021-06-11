path%% startup: add repo to path
disp(pwd)

%recursively adds all repo files/folders
addpath( genpath(pwd()) ); 

%path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
%addpath(genpath(path2repo)); %recursively adds all repo files/folders