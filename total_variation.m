function tvs = total_variation(V, L)
%% V = [v1, ..., vk] column vectors
%  L is the graph laplacian D-A
    [rows_l, cols_l] = size(L);
    if rows_l ~= cols_l
        error('laplacian must be square matrix');
    end

    [rows_v, ~] = size(V);
    if (rows_v ~= cols_l)
        error('V and Ls shape incompatible');
    end
    
    % V'LV(i,j) = vi'Lvj
    % we are interested in total variation of each vector. 
    % Total variation of v wrt to L is (1/2)*v'Lv. Thus we are 
    % only interested in the diagonal.
    tvs = diag(V' * L * V);
end