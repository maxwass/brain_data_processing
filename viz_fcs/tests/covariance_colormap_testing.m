close all;
clear;
%build colormap for covariances
path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo));
load("misc/correlations_colormap.mat");
load('dtseries_testing.mat','current_dtseries');
[num_rois, ~] = size(current_dtseries);
roi_idxs = (20:num_rois);


current_dtseries;
dtseries = current_dtseries(roi_idxs,:);
mean_signal = mean(dtseries,2);
dtseries_center = dtseries-mean_signal;
[covs, corrs, aw] = windowed_fcs(dtseries_center, 30, 10);
[N, ~, num_windows] = size(covs);

zero_diag = ones(N) - eye(N);

%extract upper triangular part with diagonal (they're symmetric)
elems_per_matrix = N*(N+1)/2;
all_vals = zeros(num_windows,elems_per_matrix);
all_vals_normed = zeros(num_windows,elems_per_matrix);
covs_normed = zeros(size(covs));

for j = 1:num_windows
    A = covs(:,:,j);%.*zero_diag;
    At = A.';
    m  = tril(true(size(At)));
    v  = At(m).';
    all_vals(j,:) = v;
    all_vals_normed(j,:) = v./norm(v); %dont include lower diag in norm!
    covs_normed(:,:,j) = A./norm(v);
end
all_vals = reshape(all_vals, [],1);
all_vals_normed = reshape(all_vals_normed, [],1);

%{
covs_zero_diag = zeros(size(covs));
covs_normed = zeros(size(covs));
for j = 1:num_windows
   cov = covs(:,:,j);
   covs_zero_diag(:,:,j) = covs(:,:,j).*zero_diag;
   
   covs_normed(:,:,j) = covs(:,:,j)./norm(covs(:,:,j));
end
%}



%plot histogram of datavalues

figure;
t = tiledlayout(1,2);
nexttile(t)
histogram(all_vals, 'Normalization', 'probability');
title('RAW COV')
nexttile(t)
histogram(all_vals_normed, 'Normalization', 'probability');
title('Normed (tri) COV')


% include pct of values
%val_pct = prctile(reshape(covs_normed(


figure;
ax_arr = gobjects(5,5);
tiled_axes = tiledlayout(5,5, 'TileSpacing','compact','Padding','none');
title(tiled_axes, 'RAW COV')
cmapLim = 2000;
for row = 1:5
    for col = 1:5
    ax = nexttile(tiled_axes);
    ax_arr(row, col) = ax;
    i = randsample(N,1);
    imagesc(covs(:, :, i));
    set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
    xlim(ax,[1,N]); xlim(ax,'manual');
    ylim(ax,[1,N]); ylim(ax, 'manual');
    daspect(ax,[1 1 1]);
    colormap(ax, colorMap);
    %set(gca,'CLim',[-cmapLim cmapLim])
    
    end
end
colorbar('NorthOutside')

cls = linspace(.05,.18);
%for cl = .04:.005:.08
for cl = 2059:300:3355
    set(ax_arr,'CLim',[-cl cl])
    txt = sprintf('clim: %f', cl);
    title(tiled_axes, txt)   
end


