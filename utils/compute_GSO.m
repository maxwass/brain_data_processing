function [S] = compute_GSO(A, GSO)
%% Contruct alternate GSOs from adjacency matrix
D_vec  = sum(A,2);
D = diag(D_vec);
D_norm = diag(D_vec.^(-.5));

if any(~D_vec)
    warning('This graph has %d nodes with 0 neighbors. Setting normalized row/col to 0', sum(~D_vec));
end

if isequal(GSO, 'L')
    L = D-A;
    S = L;
elseif isequal(GSO, 'L_norm')
    L = D-A;
    L_norm = D_norm*L*D_norm;
    % by defn of L_norm
    % (https://en.wikipedia.org/wiki/Laplacian_matrix#Symmetric_normalized_Laplacian),
    % this position should be 0
    L_norm(isnan(L_norm)) = 0.0;
    S = L_norm;
elseif isequal(GSO, 'A')
    S = A;
elseif isequal(GSO, 'A_norm')
    A_norm = D_norm*A*D_norm; % Should A_norm be 1/eig_max * A??
    A_norm(isnan(A_norm)) = 0.0;
    S = A_norm;
else
    error("Unrecognized GSO %s", GSO);
end

end


