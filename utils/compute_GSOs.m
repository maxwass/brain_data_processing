function [A_norm, L, L_norm] = compute_GSOs(A)
%% Contruct alternate GSOs from adjacency matrix

D_vec  = sum(A,2);
D = diag(D_vec);
D_norm = diag(D_vec.^(-.5));

L = D-A;
L_norm = D_norm*L*D_norm;

% if some node i has NO edges => Dii = 0 => 1/sqrt(0) => Inf => Nan 

% by defn of L_norm
% (https://en.wikipedia.org/wiki/Laplacian_matrix#Symmetric_normalized_Laplacian),
% this position should be 0
L_norm(isnan(L_norm)) = 0.0;

A_norm = D_norm*A*D_norm;
A_norm(isnan(A_norm)) = 0.0;

if any(~D_vec)
    warning('This graph has %d nodes with 0 neighbors. Setting normalized row/col to 0', sum(~D_vec));
end


end

