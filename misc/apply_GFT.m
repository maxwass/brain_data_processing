function [signals_freq, GFT, evals] = apply_GFT(signals,subject,atlas, include_subcortical, GSO)
%Applies GFT to matrix of signals (each column is a signal)
% subject = string subject identifier

    
%% load scs
if atlas=="desikan"
    scs_file = load('scs_desikan.mat');
elseif atlas=="destrieux"
    scs_file = load('scs_desikan.mat');
else
    error("NOT IMPLIMENTED")
end
    

%% find index of subject in tensor
idx = find(scs_file.subject_list==int64(str2num(subject)));
A_full_raw = scs_file.scs(:,:,idx);
A_full_transform = log(A_full_raw+A_full_raw' +1); %transform reccomended by Zhengwu

A_cortical_transform = A_full_transform(20:end, 20:end);

if include_subcortical == 1
    A = A_full_transform;
else
    A = A_cortical_transform;
end

%% compute GFT and freqs
D_vec  = sum(A,2);
D = diag(D_vec);
D_norm = diag(D_vec.^(-.5));
L = D-A;
L_norm = D_norm*L*D_norm;
A_norm = D_norm*A*D_norm;

if strcmp(GSO, 'A')
    [evecs, evals] = eig(A);
elseif strcmp(GSO, 'A_norm')
    [evecs, evals] = eig(A_norm);
elseif strcmp(GSO, 'L')
    [evecs, evals] = eig(L);
elseif strcmp(GSO, 'L_norm')
    [evecs, evals] = eig(L_norm);
else
    error('unrecognized GSO requested!%s', GSO)
end

%eigenvalues in sorted order
[d,ind] = sort(diag(evals));
evals = evals(ind,ind);
evecs = evecs(:,ind);

GFT = transpose(evecs);


%% apply GFT
signals_freq = GFT*signals;

clear scs_file
end

