function plot_scalars(ax, signals, signals_freq, covs, corrs, S)
%% Compute and plot scalars to plot on single axes

%% 
[N, ~, num_windows] = size(covs);
[num_evals, num_obsvs] = size(signals_freq);

%% define scalar metrics to use for plotting
cov_energy           = zeros(num_windows, 1);
cov_energy_wo_diag   = zeros(num_windows, 1);
inner_product    = zeros(num_windows, 1);
signals_energy = zeros(num_windows, 1);


high_freq_cutoffs = 10;
low_freq_cutoffs  = 34;
med_freq_cutoffs  = 10;
high_freq_ratios  = zeros(num_windows,high_freq_cutoffs);
low_freq_ratios   = zeros(num_windows,low_freq_cutoffs);
med_freq_ratios   = zeros(num_windows,med_freq_cutoffs);


%% loop through, calculate each quantity, and store...
for j = 1:num_windows
    x    = signals(:,j);
    cov  = covs(:,:,j);
    corr = corrs(:,:,j);
    
    signals_energy(j) = norm(x,2);
    all_cov_vals         = extract_upper_tri(cov,0);
    all_cov_vals_wo_diag = extract_upper_tri(cov,1);
    
    cov_energy(j)         = norm(all_cov_vals);
    cov_energy_wo_diag(j) = norm(all_cov_vals_wo_diag);
    
    % simple matching - is there a strong corresponence between edges and
    % corrrelations?
    %inner_product(j) = sum(corr.*S, 'all');
    
    freq_repr = signals_freq(:,j);
    %compute low freq ratios
    for num_freqs_include = 1:low_freq_cutoffs
        low_freq_energy = norm( freq_repr(1:num_freqs_include) );
        energy_ratio    = low_freq_energy/ norm( freq_repr );
        low_freq_ratios(j, num_freqs_include) = energy_ratio;
    end
    
    %compute med freq ratios
    middle_idx = ceil(num_evals/2);
    for freq_idx = 0:(med_freq_cutoffs-1)
        range = (middle_idx-freq_idx):(middle_idx+freq_idx); %size 1, 3, 5, 7,...
        med_freq_energy = norm( freq_repr(range) );
        energy_ratio    = med_freq_energy/ norm( freq_repr );
        med_freq_ratios(j, freq_idx+1) = energy_ratio;
    end

    %compute high freq ratios
    for num_freqs_include = 1:high_freq_cutoffs
        start_eig = num_evals-num_freqs_include+1; %==num_evals when num_freqs_include = 1
        high_freq_energy = norm( freq_repr(start_eig:end) );
        energy_ratio = high_freq_energy/ norm( freq_repr );
        high_freq_ratios(j, num_freqs_include) = energy_ratio;
    end
    
end

%% Only plot line best correlated (pos or neg) with *energy* <- CHANGE THIS
%{
direction = 'abs';

% which # low freqs is best correlated with energy?
[best_corr_low, best_corr_low_p, best_low_idx] = most_correlated(cov_energy, low_freq_ratios, direction);
low_freq_2_plot = low_freq_ratios(:, best_low_idx);
num_freqs_low   = best_low_idx;

% which # med freqs is best correlated with energy?
[best_corr_med, best_corr_med_p, best_med_idx] = most_correlated(cov_energy, med_freq_ratios, direction);
med_freq_2_plot = med_freq_ratios(:, best_med_idx);
num_freqs_med   = (best_med_idx*2)-1;

% which # high freqs is best correlated with energy?
[best_corr_high, best_corr_high_p, best_high_idx] = most_correlated(cov_energy, med_freq_ratios, direction);
high_freq_2_plot = high_freq_ratios(:, best_high_idx);
num_freqs_high   = best_high_idx;
%}

%% perform plotting on given ax
yyaxis(ax,'left'); %main metric on left
plot(ax, (cov_energy-mean(cov_energy))/100, '-k', 'LineWidth', 1, 'DisplayName', 'windowed cov energy/100');
hold(ax, 'on');
plot(ax, (signals_energy-mean(signals_energy)), '-b', 'LineWidth', 1, 'DisplayName', 'windowed signal energy');

%plot(ax, energy_wo_diag-mean(energy_wo_diag), '--k', 'LineWidth', 1, 'DisplayName', 'energy w/o diag (self cov)');
set(ax, 'YColor', 'k') % set y axis color to same as line color
ylabel(ax, 'Energy (mean centered)');

yyaxis(ax,'right');
num_freq_low = low_freq_cutoffs;
low_freq_2_plot = low_freq_ratios(:,num_freq_low);
txt_low = sprintf('%d (of %d) lowest freqs', int16(num_freq_low), num_evals);

num_freq_med = med_freq_cutoffs;
med_freq_2_plot = med_freq_ratios(:,num_freq_med);
txt_med = sprintf('%d (of %d) middle freqs', int16(num_freq_med*2-1), num_evals);

num_freq_high = high_freq_cutoffs;
high_freq_2_plot = high_freq_ratios(:,num_freq_high);
txt_high = sprintf('%d (of %d) highest freqs', int16(num_freq_high), num_evals);

plot(ax, 100*low_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_low);
%plot(ax, 100*med_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_med);
%plot(ax, 100*high_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_high);

%{
% plot those subsets which correlate best with cov energy
txt_low = sprintf('%d lowest freqs, corr w/ energy %.2f, p-val %.2f', ...
    int16(num_freqs_low),best_corr_low, best_corr_low_p);
txt_med = sprintf('%d medium freqs, corr w/ energy %.2f, p-val %.2f', ...
    int16(num_freqs_med),best_corr_med, best_corr_med_p);
txt_high = sprintf('%d high freqs, corr w/ energy %.2f, p-val %.2f', ...
    int16(num_freqs_high),best_corr_high, best_corr_high_p);
plot(ax, 100*low_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_low);
plot(ax, 100*med_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_med);
plot(ax, 100*high_freq_2_plot,'LineWidth', 1, 'DisplayName', txt_high);
%}
hold(ax, 'off');
ylabel(ax, '% energy in freq range', 'FontSize', 15);% (||x_high freq|| / ||x||)^2'
xlabel(ax, 'Window')
ylim(ax, [0,100])
legend(ax)

end
