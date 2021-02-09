function [lim] = find_cov_clim(cov_tensor, percentile)

[N, ~, num_windows] = size(cov_tensor);

zero_diag = ones(N) - eye(N);
%extract upper triangular part with diagonal (they're symmetric)
elems_per_matrix = N*(N+1)/2;
all_vals = zeros(num_windows,elems_per_matrix);

for j = 1:num_windows
    A = cov_tensor(:,:,j);%.*zero_diag;
    At = A.';
    m  = tril(true(size(At)));
    v  = At(m).';
    all_vals(j,:) = v;
end
all_vals = reshape(all_vals, [],1);

lim = prctile(abs(all_vals), percentile);

end

