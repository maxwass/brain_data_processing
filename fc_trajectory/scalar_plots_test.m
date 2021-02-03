% scalar plot testing
close all; clear;

path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo));

%load in first patient which is saved on machine for faster testing
subject = 100206;
load('dtseries_testing.mat','current_dtseries');
atlas = "desikan";
include_subcortical = false;

[num_rois, ~] = size(current_dtseries);
roi_idxs = (20:num_rois);
if include_subcortical
    roi_idxs = 1:num_rois;
end

%% mean center raw signals
dtseries            = current_dtseries(roi_idxs,:);
mean_signal         = mean(dtseries,2);
dtseries_center     = dtseries - mean_signal;


%% window signal
windowsize = 20;
movesize = 5;

[covs, corrs, aw]   = windowed_fcs(dtseries_center, windowsize, movesize);
[N, ~, num_windows] = size(covs);

%% load scs from saved tensor file
GSO = 'L';
A = extract_sc(subject, atlas, include_subcortical);
A_upp_tri = triu(A,0);
A_upp_tri_norm = A_upp_tri./norm(A_upp_tri,2); %only consider upper triangular part


%% define scalar metrics to use for plotting
energy           = zeros(num_windows, 1);
energy_wo_diag   = zeros(num_windows, 1);
total_covariance = zeros(num_windows, 1);
inner_product    = zeros(num_windows, 1);

high_freq_cutoffs = 10;
low_freq_cutoffs = 10;
med_freq_cutoffs = 10;
high_freq_ratios = zeros(num_windows,high_freq_cutoffs);
low_freq_ratios = zeros(num_windows,low_freq_cutoffs);
med_freq_ratios = zeros(num_windows,med_freq_cutoffs);

[signals_freq, GFT, evals] = apply_GFT(dtseries_center,subject, atlas, include_subcortical, GSO);

%% loop through, calculate each quantity, and store...
zero_diag = ones(N) - eye(N);

%extract upper triangular part with diagonal (they're symmetric)
elems_per_matrix = N*(N+1)/2;
covs_normed = zeros(size(covs));

for j = 1:num_windows
    cov  = covs(:,:,j);
    corr = corrs(:,:,j);
    
    all_cov_vals         = extract_upper_tri(cov,0);
    all_cov_vals_wo_diag = extract_upper_tri(cov,1);
    
    energy(j)           = norm(all_cov_vals);
    energy_wo_diag(j)   = norm(all_cov_vals_wo_diag);
    total_covariance(j) = sum(all_cov_vals);
    
    % simple matching - is there a strong corresponence between edges and
    % corrrelations?
    inner_product(j)   = sum(corr.*A_upp_tri_norm, 'all');
    
    freq_repr = signals_freq(:,j);
    %compute high freq ratios
    for num_freqs_include = 1:high_freq_cutoffs
        start_eig = length(evals)-num_freqs_include+1; %==length(evals) when num_freqs_include = 1
        high_freq_energy = norm( freq_repr(start_eig:end) );
        energy_ratio = high_freq_energy/ norm( freq_repr );
        high_freq_ratios(j, num_freqs_include) = energy_ratio;
    end
    
    %compute low freq ratios
    for num_freqs_include = 1:low_freq_cutoffs
        low_freq_energy = norm( freq_repr(1:num_freqs_include) );
        energy_ratio    = low_freq_energy/ norm( freq_repr );
        low_freq_ratios(j, num_freqs_include) = energy_ratio;
    end
    
    %compute med freq ratios
    middle_idx = ceil(length(evals)/2);
    for freq_idx = 0:(low_freq_cutoffs-1)
        range = (middle_idx-freq_idx):(middle_idx+freq_idx); %1, 3, 5, 7,...
        med_freq_energy = norm( freq_repr(range) );
        energy_ratio    = med_freq_energy/ norm( freq_repr );
        med_freq_ratios(j, freq_idx+1) = energy_ratio;
    end
end

%correlations

%% which # low freqs is best correlated with energy?
[e_low_corrs, e_low_ps] = corrcoef([energy,low_freq_ratios]);

%first row (2:end) has energy correlations with low_freq_ratios
a = e_low_corrs(1,:);
b = a(2:end);
[best_corr_low, best_corr_idx] = max(abs(b));
best_corr_low_p_val = e_low_ps(1,(best_corr_idx+1));
low_freq_2_plot = low_freq_ratios(:, best_corr_idx);
num_freqs_low = best_corr_idx;

%% which # med freqs is best correlated with energy?
[e_med_corrs, e_med_ps] = corrcoef([energy,med_freq_ratios]);

%first row (2:end) has energy correlations with low_freq_ratios
a = e_med_corrs(1,:);
b = a(2:end);
[best_corr_med, best_corr_idx] = max(abs(b));
best_corr_med_p_val = e_med_ps(1,(best_corr_idx+1));
med_freq_2_plot = med_freq_ratios(:, best_corr_idx);
num_freqs_med = (best_corr_idx-1)*2;

%% which # high freqs is best correlated with energy?
[e_high_corrs, e_high_ps] = corrcoef([energy,high_freq_ratios]);

%first row (2:end) has energy correlations with low_freq_ratios
a = e_high_corrs(1,:);
b = a(2:end);
[best_corr_high, best_corr_idx] = max(abs(b));
best_corr_high_p_val = e_high_ps(1,(best_corr_idx+1));
high_freq_2_plot = high_freq_ratios(:, best_corr_idx);
num_freqs_high = best_corr_idx;






figure;
ax = axes;
yyaxis(ax,'left');
h1 = plot(energy-mean(energy), '-k', 'LineWidth', 1, 'DisplayName', 'energy'); 
set(gca, 'YColor', 'k') % set y axis color to same as line color
hold on;
plot(energy_wo_diag-mean(energy_wo_diag), '--k', 'LineWidth', 1, 'DisplayName', 'energy w/o diag (self cov)'); hold on;
%plot(x, total_covariance, ':k', 'LineWidth', 1, 'DisplayName', 'sum of cov entries'); hold on;
hold off;
ylabel('Energy');

yyaxis(ax,'right');
%plot(x,inner_product, 'b', 'LineWidth', 1, 'DisplayName', 'inner product: <A,cov>'); 
%{
for es = 1:high_freq_cutoffs
    fr = 100*high_freq_ratios(:,es);
    plot(100*high_freq_ratios(:,es),'LineWidth', 1, 'DisplayName', ...
        sprintf('%d largest freqs', int16(es)));
    hold on;
    
end
%}
txt_low = sprintf('%d lowest freqs, corr w/ energy %.2f, p-val %.2f', ...
    int16(num_freqs_low),best_corr_low, best_corr_low_p_val);
txt_med = sprintf('%d medium freqs, corr w/ energy %.2f, p-val %.2f', ...
    int16(num_freqs_med),best_corr_med, best_corr_med_p_val);
txt_high = sprintf('%d medium freqs, corr w/ energy %.2f, p-val %.2f', ...
    int16(num_freqs_high),best_corr_high, best_corr_high_p_val);
plot(100*low_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_low); hold on;
plot(100*med_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_med); hold on;
plot(100*high_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_high); hold off;
ylabel('% energy', 'FontSize', 15);% (||x_high freq|| / ||x||)^2'
xlabel('Window')
ylim([0,100])
legend

%scalar_observations = [energy; energy_wo_diag; high_freq_ratios];
%scalar_corrs = corrcoef(scalar_observations);