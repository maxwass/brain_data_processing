function [GFT, evals_vec, S] = extract_GFT(subject, atlas, include_subcortical, GSO)
%Construct GFT from GSO of subject's sc scan
% subject = string subject identifier

%load scs from saved tensor file
A = extract_sc(subject, atlas, include_subcortical);

%% compute GFT and freqs
[S] = compute_GSO(A, GSO);
[evecs, evals] = eig(S);

%eigenvalues in sorted order
[~,ind] = sort(diag(evals));
evals = evals(ind,ind);
evals_vec = diag(evals);
evecs = evecs(:,ind);

%ensure eigenvector corresponnding to 0 freq (constant) is positive. For
%laplacian, this represents the average value of all nodes, which for us
%is non-negative. Thus for more intuitive plotting, make positive.
zero_idx = find( abs(evals_vec) < 1e-8 );
if ~isempty(zero_idx)
    
    % this should likely be equal to zero (Laplacian, Laplacian norm)
    evals_vec(zero_idx) = 0.0;
    
    % pointing in 'positive' (1st quadrant) or 'negative' (3rd quadrant) direction?
    zero_eigvec = evecs(:, zero_idx);
    is_neg      = dot(ones(size(zero_eigvec)), zero_eigvec) < 0;
    
    % if eigvector pointing in negative direction, make positive
    if is_neg
        evecs(:, zero_idx) = -1*evecs(:, zero_idx);
    end
end

    
GFT = transpose(evecs); % GFT == V^H
end

