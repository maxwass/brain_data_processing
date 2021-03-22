function zcs = zero_crossings(V, A)
%% V is a matrix of signals (column vector) in R^n supported on graph G
%  A is the adjacency matrix in R^(n x n) of graph G

    [rows_a, cols_a] = size(A);
    if rows_a ~= cols_a
        error('adjacency must be square matrix');
    end

    [rows_v, ~] = size(V);
    if (rows_v ~= cols_a)
        error('V and Ls shape incompatible');
    end
    
    zcs = zero_crossings_from_total_variation(V,A);
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
    
    %if undirected graph, divide by two
    if abs(A-A')<.000001
        zcs = zcs/2;
    end
end

function zcs = zero_crossings_from_total_variation(V,A)
    zcs = (1/4) * total_variation(sign(V),sign(A));
end

function zcs_test()
% Copy and Paste this into the command line
A = [[0,1,0];[1,0,1];[0,1,0]];
z1 = [1,-1,-1]';    soln1 = 1;
z2 = [2,1,2]';      soln2 = 0;
z3 = [0,0,0]';      soln3 = 0;
z4 = [3.3,-2.1,1]'; soln4 = 2;


op = @(z,A) zero_crossings_from_outer_product(z,A);
tv = @(z,A) zero_crossings_from_total_variation(z,A);
assert(op(z1,A)==tv(z1,A)==soln1);
assert(op(z2,A)==tv(z2,A)==soln2);
assert(op(z3,A)==tv(z3,A)==soln3);
assert(op(z4,A)==tv(z4,A)==soln4);



end
