function [signals_freq, GFT, evals] = apply_GFT(signals,subject,atlas, include_subcortical, GSO)
%Applies GFT to matrix of signals (each column is a signal)
% subject = string subject identifier


%load scs from saved tensor file
A = extract_sc(subject, atlas, include_subcortical);

%% compute GFT and freqs
[A_norm, L, L_norm] = compute_GSOs(A);

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
[~,ind] = sort(diag(evals));
evals = evals(ind,ind);
evecs = evecs(:,ind);

GFT = transpose(evecs); % GFT == V^H


%% apply GFT
signals_freq = GFT*signals;
end
