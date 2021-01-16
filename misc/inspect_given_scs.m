
%scs
% subject_list (1065x1 int64)
% scs (87x87x1065 double)
sc_file = load('~/Documents/MATLAB/brain_data_preprocess/data/scs_desikan.mat');
SCs = sc_file.scs; 
subject_list_sc = sc_file.subject_list;

%any graphs nan?
anynan = any(any(any(isnan(SCs))));

%any matrix have all zeros?
allzero = all(any(any(SCs>0,2)));

%any negative?
anyneg = any(any(any(SCs<0,2)));

if anynan || allzero || anyneg
    fprintf('bad SC graph!')
end

%viz
num_scs = length(SCs);
for i = 1:num_scs
    disp(i)
    disp(SCs(1:3,1:3,i))
end