%% How does sparsity change if we remove small/large edges?

clear; close;
atlas = 'desikan';
include_subcortical = false;

scs_upper_tri = load_processed_scs(atlas, include_subcortical, 'log');
scs = apply_to_tensor_slices(@(x) (x+x'), scs_upper_tri);


thresholds = 0:.1:10;
mean_sparsities = zeros(length(thresholds), 1);
median_sparsities = zeros(length(thresholds), 1);

for idx = 1:length(thresholds)
    sparsities  = compute_sparsities(scs_upper_tri, thresholds(idx), 10);
    mean_sparsities(idx) = mean(sparsities);
    median_sparsities(idx) = median(sparsities);  
end


t = tiledlayout(1,1);
ax = nexttile();
plot(ax, thresholds, mean_sparsities, 'DisplayName', 'Mean Sparsity');
hold on;
plot(ax, thresholds, median_sparsities, 'DisplayName', 'Mean Sparsity');

xlabel('Threshold', 'FontSize', 20);
ylabel('Sparsity', 'FontSize', 20);
title('Sparsity vs Threshold', 'FontSize', 35);
legend;


% input must be 3D tensor: N x N x num_matrices
function sparsities  = compute_sparsities(scs_upper_tri, low_threshold, high_threshold)
    [N, ~, ~] = size(scs_upper_tri);
    num_edges_fully_connected = N*(N-1)/2; % no self loop, undirected

    % perform sum on each slice.
    num_edges = squeeze(sum(low_threshold<scs_upper_tri & scs_upper_tri<high_threshold, [1,2]));
    
    sparsities = num_edges/num_edges_fully_connected;
end

function scs_upper_tri = load_processed_scs(atlas, include_subcortical, edge_transform)
    cortical_idxs = get_roi_idxs(atlas, include_subcortical);

    if isequal(atlas, 'desikan')
        scs_file = load('scs_desikan.mat');
    else
        error('atlas %s unsupported', atlas)
    end
    
    scs_upper_tri = scs_file.scs(cortical_idxs, cortical_idxs, :);

    if isequal(edge_transform, 'log')
        scs_upper_tri = log(scs_upper_tri+1);
    end


end

