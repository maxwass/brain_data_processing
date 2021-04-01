%% sparsity distribution of SCs
clear; close;
atlas = 'desikan';
include_subcortical = false;

scs_upper_tri = load_processed_scs(atlas, include_subcortical, 'log');
scs = apply_to_tensor_slices(@(x) (x+x'), scs_upper_tri);

sparsities  = compute_sparsities(scs_upper_tri, 0);

t = tiledlayout(1,1);
ax = nexttile(t);
histogram(ax, sparsities, 'DisplayName', 'SC sparsities');
%xlim(ax, [0,1]);
%xline(ax, 0.5, '--b', 'LineWidth', 3);
ylabel('patient counts', 'FontSize', 20);

hold on;
yyaxis right;
[~, stats] = cdfplot(sparsities);
ylabel('Pr(s < sparsity)', 'FontSize', 20);
xlabel('s', 'FontSize', 20);
stats_text = sprintf("min: %.2f \nmax: %.2f \nmean: %.2f \nmedian: %.2f \nstd: %.2f", stats.min, stats.max, stats.mean, stats.median, stats.std);
text(.2, .6, stats_text, 'Units','normalized', 'FontSize',25);

if ~include_subcortical
    txt = 'w/ only Cortical Nodes';
else
    txt = 'w/ Cortical and Subcortical Nodes';
end
title_txt = sprintf('SC sparsities. %s Atlas %s', atlas, txt);
title('SC sparsities. Desikan Atlas w/ only Cortical Nodes.', 'FontSize', 30);
hold off;

% sort SCs by sparsities
[~, I] = sort(sparsities, 'ascend');
scs_upper_tri = scs_upper_tri(:, :, I);
scs = scs(:, :, I);

% input must be 3D tensor: N x N x num_matrices
function sparsities  = compute_sparsities(scs_upper_tri, threshold)
    [N, ~, ~] = size(scs_upper_tri);
    num_edges_fully_connected = N*(N-1)/2; % no self loop, undirected

    % perform sum on each slice.
    num_edges = squeeze(sum(scs_upper_tri>threshold, [1,2]));
    
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
