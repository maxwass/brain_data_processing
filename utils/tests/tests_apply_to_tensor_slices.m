%% testing scrupt for get_slice_idxs and apply_to_tensor_slices

% 
%% get_slice_idxs
% give tensor get indices for ith 'slice' of that data
% input: 2x5 -> ith slice is column i -> {':', [i]}
two_d = @(i) {':', i};
three_d = @(i) {':', ':', i};
% input: 3x3x5 -> ith slice is submatrix i -> {':',':',[i]}
% input: 3x1x3 -> ith slice is column/matrix i -> {':',':',[i]}
% input: 1x3x3 -> ith slice is row/matrix i -> {':',':',[i]}

sz1 = size(ones(2,5));
sz2 = size(ones(3,3,5));
sz3 = size(ones(3,1,3));
sz4 = size(ones(1,3,3));

inputs = {sz1,sz2,sz3,sz4};
truths = {two_d,three_d,three_d,three_d};

for i = 1:length(inputs)
    in_size = inputs{i};
    truth_out = truths{i};
    for j = 1:in_size(end) %check all slices
        tensor1 = get_slice_idxs(length(in_size), j);
        b = truth_out(j);
        for k = 1:length(tensor1)
            if tensor1{k}~=b{k}
                error('%dth input is inccorect', i);
            end
        end
    end
end


%% apply_to_tensor_slices

% apply function to each slice
tensor1 = .5.*ones(3,3,2);
tensor2 = ones(3,5,2);
tensor3 = .25*ones(1,3,2);
tensor4 = 2*ones(3,1,2);
ts = {tensor1, tensor2, tensor3, tensor4};

f1 = @(x) x.^2;
f2 = @(x) sqrt(x);
% psuedo-function partial application through anonymous functions
% https://stackoverflow.com/questions/9154271/partial-function-evaluation-in-matlab
norm_spec = @(spec,x) @(x) norm(x,spec);
funcs = {f1, f2, norm_spec(2), norm_spec(Inf), norm_spec('fro')};

% MUST MANUALLY CHECK
for i = 1:length(funcs)
    f = funcs{i};
    atts = @(x) apply_to_tensor_slices(f,x);
    outs = cellfun(atts, ts, 'UniformOutput',false);
    
end

a_out = apply_to_tensor_slices(norm_spec('fro'), tensor3);

data = rand(3,100);
c = cov(data');
c = repmat(c,1,1,2);
a_out = apply_to_tensor_slices(@corrcov, c);

