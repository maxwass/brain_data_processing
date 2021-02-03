function [v] = extract_upper_tri(A,k)
%return all elements on and above the kth diagonal of A as a vector
m  = triu(true(size(A)),k);
v  = A(m);
end