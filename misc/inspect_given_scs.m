
load('HCP_subcortical_CMData_desikan.mat')
SCs = squeeze(loaded_tensor_sub(:,:,1,:)); % (87x87x1065 double)

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