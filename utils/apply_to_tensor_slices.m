function [out_tensor] = apply_to_tensor_slices(f, in_tensor)
%% tensor :: N x K x M
%  f is a function that takes a slice of in_tensor and outputs a new
%  slice/scalar
%   - to call f use normally f(input)
%   - when providing a function as input use @:
%   apply_to_tensor_slices(@norm, A);
%% f: N x K array -> _ array
%  ex) matrix multiply signals
%     -tensor :: 68 x 68 x 120 array. f(tensor_slice) :: 68 x K array. output :: 68 x K x 120 array
%% if tensor: N x M
%  f: N x M -> _ array (apply f to columns)
%  ex) compute energies
%     -tensor :: 68 x 120 array. f(tensor_slice) :: 1 x 1 array. output :: 1 x 120 array
%  ex) compute freq repr
%     -tensor :: 68 x 120 array. f(tensor_slice) :: 68 x 1. output :: 68  x 120 array

% for information on variable indexing using cell arrays:
% https://www.mathworks.com/matlabcentral/answers/362211-variable-indexing-for-n-dimension-data


% input tensor info
sz = size(in_tensor);
M = sz(end); %number of slices

% infer output size from output size of f?
sample_in  = ones(sz(1:end-1));
sample_out = f(sample_in);
out_size = cat(2, size(sample_out), M);

% init output tensor of size out_size_f x M
out_tensor = zeros(out_size);

for which_slice = 1:M
    in_slice_idxs = get_slice_idxs(length(size(in_tensor)), which_slice);
    in_slice      = in_tensor(in_slice_idxs{:});    % in_tensor((:,:,2))
    
    
    out_slice_idxs = get_slice_idxs(length(size(out_tensor)), which_slice);
    out_tensor(out_slice_idxs{:}) = f(in_slice);
end

%squeeze out_tensor to get rid of possible singleton dimension??

