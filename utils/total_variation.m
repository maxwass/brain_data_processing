function tvs = total_variation(V, A)
%% V = [v1, ..., vk] column vectors
%  A is the graph adjacency
    [rows_a, cols_a] = size(A);
    if rows_a ~= cols_a
        error('adjacency must be square matrix');
    end

    [rows_v, ~] = size(V);
    if (rows_v ~= cols_a)
        error('V and Ls shape incompatible');
    end
    %  L is the graph laplacian D-A
    L = diag(sum(A,2)) - A;
    
    % V'LV(i,j) = vi'Lvj
    % we are interested in total variation of each vector. 
    % Total variation of v wrt to L is (1/2)*v'Lv. Thus we are 
    % only interested in the diagonal.
    tvs = diag(V' * L * V);
end

function test()
    % Copy and Paste this into the command line

    A1 = [[0,1,0];[1,0,1];[0,1,0]];
    z1 = [1,1,1]; soln1 = 0;
    z2 = [0,0,0]; soln2 = 0;
    z3 = [2.5,-1.5,3.5]; soln3 = 41;
    V = [z1;z2;z3]';
    tvs = total_variation(V,A1);
    assert(all(tvs==[soln1,soln2,soln3]'));
    
    A2 = A1;
    A2(1,2) = 2;
    A2(2,1) = 2;
    soln3 = 57;
    tvs = total_variation(V,A2);
    assert(all(tvs==[soln1,soln2,soln3]'));
    
end