function [A_norm, L, L_norm] = compute_GSOs(A)
%% Contruct alternate GSOs from adjacency matrix

D_vec  = sum(A,2);
D = diag(D_vec);
D_norm = diag(D_vec.^(-.5));

L = D-A;
L_norm = D_norm*L*D_norm;

A_norm = D_norm*A*D_norm;


end

