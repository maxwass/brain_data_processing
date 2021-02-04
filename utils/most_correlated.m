function [best_corr, best_corr_p, best_idx] = most_correlated(x, y, direction)
%% given data vector x in Mx1 and data matrix y in MxK, compute pearson correlation
% between x and each column of y. Find 'best' correlation (positive,
% negtive, or abs linear relationship), and return the correlation, p-val,
% and column of y which corresponds.
if length(size(x))~=2 || (length(size(y)) ~= 2)
    error('x and/or y has wrong input size. x is %d-dim, y is %d-dim\n', length(size(x)), length(size(y)));
end

[Mx, a] = size(x);
[My, ~] = size(y);
if Mx~=My || a~=1
    error('x and y have incompatible sizes.');
end


%% which # high freqs is best correlated with energy?
[cs, ps] = corrcoef([x,y]);

%extract first row, ignoring self correlation (first element)
x_cs = cs(1,2:end);

if strcmp(direction, 'neg')
    [best_corr, best_idx] = min(x_cs); %take most negtive correlation
elseif strcmp(direction, 'abs')
    [best_corr, best_idx] = max(abs(x_cs)); %take strongest linear relationship (pos or neg)
else
    [best_corr, best_idx] = max(x_cs); %take most positive correlation
end

best_corr_p = ps(1,(best_idx+1)); %add one bc removed first element first row of cs

end

