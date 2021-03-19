function zcs = zero_crossings(z, A, impl)
%% z is a column vector in R^n supported on graph G
%  A is the adjacency matrix in R^(n x n) of graph G
%  impl is a string specifiying which implimentation to use

    [rows_a,cols_a] = size(A);
    if rows_a ~= cols_a
        error('A must be square: %d x %d', rows_a ,cols_a);
    end
    
    % z must be vector
    [rows_z, cols_z] = size(z);
    if ~(rows_z==1 || cols_z==1)
        error('z must be a vector: %d x %d', rows_z, cols_z);
    elseif cols_z > rows_z
        % given row vector, use column vector
        z = z';
    end
    
    % dimensions of z and A must match
    if rows_z ~= cols_a
        error('z (%d, %d) and A (%d, %d) do not have compatible shape', rows_z, cols_z, rows_a, cols_a);
    end
    
    %% which implimentation should we use?
    if isequal(impl, 'outer_product')
        zcs = (1/2) * zero_crossings_from_outer_product(z,A);
    else
        zcs = zero_crossings_from_total_variation(z,A);
    end
    
    %zcs_test()
end

function zcs = zero_crossings_from_outer_product(z,A)
% z*z' is an R^nxn matrix.
%   -the ith row of z*z' is the scalar product of z(i) with all other entries
%     (including itself).
%   - negative entries in the ith row, say for example z(i,j)=-.35, indicate that z(i) and z(j) have differnt sign.
%   - z*z'< 0 sets all such negative entries to 1, everything else set to 0.
% A is the adjacency matrix supporting signal z.
%   - A>0 simply tells us which nodes have an edge.
% We are only interested in differing signs between signal elements
%   connected by an edge.
% (z*z'<0) .* (A>0) 0's elements NOT connected by edge. Doesn't affect others.
%   -the number of nonzeros entries in this matrix is the number of zero
%   crossings. 
% For undirected graphs we can divide by 2. Symmetric relationship.
    zcs = sum( (z*z'<0).*(A>0), 'all');
end

function zcs = zero_crossings_from_total_variation(z,A)
    A_sign = sign(A); 
    L_sign = diag(sum(A_sign,2)) - A_sign;
    zcs = (1/4) * total_variation(sign(z),L_sign);
end

function zcs_test()

A = [[0,1,0];[1,0,1];[0,1,0]];
z1 = [1,-1,-1]';    soln1 = 1;
z2 = [2,1,2]';      soln2 = 0;
z3 = [0,0,0]';      soln3 = 0;
z4 = [3.3,-2.1,1]'; soln4 = 2;
% mult by .5 bc undirected graph
assert(.5*zero_crossings_from_outer_product(z1,A)==soln1);
assert(.5*zero_crossings_from_outer_product(z2,A)==soln2);
assert(.5*zero_crossings_from_outer_product(z3,A)==soln3);
assert(.5*zero_crossings_from_outer_product(z4,A)==soln4);

assert(zero_crossings_from_total_variation(z1,A)==soln1);
assert(zero_crossings_from_total_variation(z2,A)==soln2);
assert(zero_crossings_from_total_variation(z3,A)==soln3);
assert(zero_crossings_from_total_variation(z4,A)==soln4);


end
