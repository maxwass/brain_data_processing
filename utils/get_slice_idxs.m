function [slice_idxs] = get_slice_idxs(num_tensor_dims, which_slice)
%Return cell array of indices to get a slice of a tensor
% ex) given 2D tensor size N x K, return {':', which_slice} <=> column which_slice
% ex) given 3D tensor size N x K x M return {':', ':', which_slice} <=>
% matrix which_slice
    slice_idxs = cell(1,num_tensor_dims);
    slice_idxs(1,:)   = {':'};
    slice_idxs(1,end) = {which_slice};
end

