path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo)); %recursively adds all repo files/folders

%scs
% subject_list (1065x1 int64)
% scs (87x87x1065 double)
sc_file = load('data/scs_desikan.mat');
SCs = sc_file.scs; 
subject_list_sc = sc_file.subject_list;

%any graphs nan?
anynan = any(any(any(isnan(SCs))));

%any matrix have all zeros?
allzero = all(any(any(SCs>0,2)));

%any negative?
anyneg = any(any(any(SCs<0,2)));

%any nodes with no neighbors?
scs_t = apply_to_tensor_slices( @(x) log(x+x'+1), SCs);

degrees = sum(scs_t, 2);
zero_degree_slice = @(d, eps) any( abs(d-0.0)<eps );
scs_with_zero_degree = zero_degree_slice(degrees, .00001);
zero_degree_idxs = find(scs_with_zero_degree);
patients_with_zero_degree = subject_list_sc(zero_degree_idxs);
fprintf('These patients have at least one node with zero degree');
disp(patients_with_zero_degree);

%find the nodes which have 0 degree
for l = 1:length(zero_degree_idxs)
    idx = zero_degree_idxs(l);
    ds = degrees(:,:,idx);
    zero_nodes = find( abs(ds)<.00001);
    disp(zero_nodes);
end

% it is always node 8 (a subcortical node)
% lets find the median of row 8
scs_t_median = median(scs_t, 3);
row_8 = scs_t_median(scs_t_median(8,:));

    

anyzerodegree = any(zero_degree_slice(degrees, .00001));

if anynan || allzero || anyneg || anyzerodegree
    fprintf('bad SC graph!')
end

%viz
num_scs = length(SCs);
for i = 1:num_scs
    disp(i)
    disp(SCs(1:3,1:3,i))
end