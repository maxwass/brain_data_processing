%% Explore properties of SCs
clear; close;

%% load data
atlas = 'desikan';
include_subcortical = false;

scs_upper_tri_raw = load_processed_scs(atlas, include_subcortical, 'raw');
scs_upper_tri_log = load_processed_scs(atlas, include_subcortical, 'log');

%% Plot Edge weight distribution
fs = 25;
fig = figure();
t = tiledlayout(2,1);

ax = nexttile(t);
title_txt = 'Raw edge weights across all SCs';
plot_ew_hist(ax, scs_upper_tri_raw(scs_upper_tri_raw>0), 800, 25, title_txt, 1500);

ax = nexttile(t);
title_txt = 'Log(edge weight + 1) across all SCs';
plot_ew_hist(ax, scs_upper_tri_log(scs_upper_tri_log>0), 250, 25, title_txt, 10);






%{
%% Do SCs with different sparsities have different edge weight distributions?
% sort SCs by sparsities
sparsities  = compute_sparsities(scs_upper_tri, 0);
[~, I] = sort(sparsities, 'ascend');
scs_upper_tri = scs_upper_tri(:, :, I);
scs = scs(:, :, I);

group_1 = 1:round((length(scs)/10));
group_2 = round(.25*length(scs)):round(.75*length(scs));
group_3 = round((9*length(scs)/10)):length(scs);

group_1_scs = scs_upper_tri(:,:,group_1);
group_1_edges = group_1_scs(group_1_scs>0);
group_2_scs = scs_upper_tri(:,:,group_2);
group_2_edges = group_2_scs(group_2_scs>0);
group_3_scs = scs_upper_tri(:,:,group_3);
group_3_edges = group_3_scs(group_3_scs>0);


ax = nexttile();
fa = .5;
ea = .5;
histogram(ax, group_1_edges, bin_edges, 'FaceAlpha', fa, 'EdgeAlpha', ea, ...
    'FaceColor', 'black', ...
    'Normalization', normalization,...
    'DisplayName', 'Lowest 10% by sparsity');
hold on;
histogram(ax, group_2_edges, bin_edges, 'FaceAlpha', fa, 'EdgeAlpha', ea, ...
    'FaceColor', 'cyan', ...
    'Normalization', normalization,...
    'DisplayName', 'Middle 33% by sparsity');
histogram(ax, group_3_edges, bin_edges, 'FaceAlpha', fa-.1, 'EdgeAlpha', ea, ...
    'FaceColor', 'red', ...
    'Normalization', normalization, ...
    'DisplayName', 'Highest 10% by sparsity');
xlim(ax, [0, xlim_high]);
legend;
hold off;
%}



function plot_ew_hist(ax, edges, num_bins, fs, title_txt, xlim_high)
    bin_edges = linspace(-1, max(edges)+1, num_bins);
    normalization = 'count';
    histogram(ax, edges, bin_edges, ...
        'Normalization', normalization, ...
        'DisplayName', 'edge weights');
    ylabel('edge counts', 'FontSize', fs+10)
    xlim(ax, [0, xlim_high]);


    hold on;
    yyaxis right;
    [~, stats] = cdfplot(edges);
    stats_text = sprintf("min: %.2f \nmax: %.2f \nmean: %.2f \nmedian: %.2f \nstd: %.2f", stats.min, stats.max, stats.mean, stats.median, stats.std);
    text(.7, .6, stats_text, 'Units','normalized', 'FontSize', fs);
    xlabel('e', 'FontSize', fs);
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15,'FontWeight','bold')
    ylabel('CDF: Pr(e < edge weight)', 'FontSize', fs);
    title(title_txt, 'FontSize', fs+5);
    hold off;

end

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


