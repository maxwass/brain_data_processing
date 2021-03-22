%% Polynomial over SCs

clear; close;
num_patients_plot = 3;
percentile = 90; 

%% some example polynomials
og_p = [.5, .5, .2]; og_p_str = '[0.5, 0.5, 0.2]';
regress_p = (1/350)*[-.18, 2.446, -.0211]; regress_p_str = 'regress';

p2 = [0.7385, -0.0107, 3.8339e-05];  p2_str = 'low pass';% [.73, -.01, 3.8e-05]'; % low pass 1
%p3 = [1, -.0005, -.00005];  p3_str = '[1, -.0005, -.00005]'; % low pass 2
p3 = [.49, .0163, -1.3061e-04]; p3_str = 'mid pass';% [.49, .0163, -.000013]'; % mid pass
p4 = (1/2)*[0, -.005, .0001];  p4_str = 'high pass';% [0, -.005, .0001]'; % high pass
%p3 = [1, -.0001, +.00008]; p3_str = '[0, .0001, -.00008]';
coeffs_list = {regress_p, p2, p3, p4};
coeffs_list_str = {regress_p_str, p2_str, p3_str, p4_str};


t = tiledlayout(num_patients_plot, 2+length(coeffs_list)+2); %2xTV/ZC + SYNTH COVS + LR COV + RL COV


atlas = 'desikan';
include_subcortical = 'false';
GSO = 'A_norm';



%% Load data
load("correlations_colormap.mat");

data = load('FreqIntervalsKept1-87', 'data').data;
rng(10);
patient_idxs = randi([1,length(data)], 1,3 );

%% Sample white signals
% Find E[GFT*w], where w is white signal
if include_subcortical
    N = 68;
else
    N = 87;
end
mu = zeros(N,1)/sqrt(N);
sigma = eye(N);
num_samples = 1000;
random_signals = mvnrnd(mu,sigma,num_samples)';


%% Plotting
for p = 1:num_patients_plot
    
    %% find eigevals (x) and freq_mean (y)
    patient_data = data(patient_idxs(p));
    [A] = extract_sc(patient_data.subject_id, atlas, include_subcortical);
    [GFT, eigvals, S] = extract_GFT(patient_data.subject_id, atlas, include_subcortical, GSO);
    % this is the mean signal in the GFT domain from random signals
    freq_mean = mean(GFT*random_signals,2);
    
    %% find polynomials to plot over freq_mean
    x_polys = min(eigvals):max(eigvals);
    y_polys = cellfun( @(cs) polyval(flip(cs), x_polys), coeffs_list, 'UniformOutput', false);
    
    
    %% interpretation of freq components: zero crossings and total variation
    zcs = zeros(length(GFT),1); % ith position is # 0-crossings of ith freq component
    V = GFT';
    L = diag(sum(A,2)) - A;
    zcs = zero_crossings(V, A);
    tvs = total_variation(V, L);
    
    % tv/zcs vs eigenvals
    ax = nexttile([1,2]);
    yyaxis(ax,'left');
    hold on;
    plot(ax, eigvals, tvs, 'DisplayName', 'total variation');
    plot(ax, eigvals, zcs, 'DisplayName', 'zero crossing');
    xticks(ax, linspace(min(eigvals), max(eigvals), 5));
    ylabel('Measure of Variation');
    xlabel('Freqs'); 
    hold off;
    
    %% Polynomial filter used
    yyaxis(ax,'right');
    hold on;
    for ps = 1:length(y_polys)
        plot(ax, x_polys, abs(y_polys{ps}), 'DisplayName', coeffs_list_str{ps});
    end
    ylabel('Filter Magnitude');
    hold off;
    legend();
    
    %% plot covs of each polynomial filter
    %H_list = create_matrix_polynomials(coeffs_list, S);
    H_list = cellfun( @(coeffs) matrix_polynomial(coeffs, S), coeffs_list, 'UniformOutput', false);
    covs = cellfun( @(H) cov((H*random_signals)' ), H_list, 'UniformOutput', false);
    
    fc_axes = gobjects(length(covs)+2,1);
    for cs = 1:length(covs)
        ax = nexttile();
        fc_axes(cs) = ax;
        imagesc(covs{cs});
        xlim(ax,[1,length(eigvals)]); xlim(ax,'manual');
        ylim(ax,[1,length(eigvals)]); ylim(ax, 'manual');
        daspect(ax,[1 1 1]);
        colormap(ax, colorMap);
        [cl] = find_cov_clim(repmat(covs{cs},1,1,2), percentile);
        set(ax,'CLim',[-cl cl]);
        cb = colorbar(ax,'EastOutside');
        title(ax, sprintf('Cov: %s', coeffs_list_str{cs}));
    end
    
    
    %% Plot actual LR/RL covs
    lr_ax = nexttile();
    fc_axes(end-1) = lr_ax;
    imagesc(patient_data.lr_fcs(20:end,20:end))
    xlim(lr_ax,[1,length(eigvals)]); xlim(ax,'manual');
    ylim(lr_ax,[1,length(eigvals)]); ylim(ax, 'manual');
    daspect(lr_ax,[1 1 1]);
    colormap(lr_ax, colorMap)
    title(lr_ax, 'Cov LR');
    
    rl_ax = nexttile();
    fc_axes(end) = rl_ax;
    imagesc(patient_data.rl_fcs(20:end,20:end))
    xlim(rl_ax,[1,length(eigvals)]); xlim(rl_ax,'manual');
    ylim(rl_ax,[1,length(eigvals)]); ylim(ax, 'manual');
    daspect(rl_ax,[1 1 1]);
    colormap(rl_ax, colorMap)
    title(rl_ax, 'Cov LR');
    
    fc_tensor = zeros(N,N,length(covs)+2);
    for fc_i = 1:length(covs)
        fc_tensor(:,:,fc_i) = covs{fc_i};
    end
    fc_tensor(:,:,end-1) = patient_data.lr_fcs(20:end,20:end);
    fc_tensor(:,:,end) = patient_data.rl_fcs(20:end,20:end);
   
    [cl] = find_cov_clim(fc_tensor(:,:,(end-1):end), percentile);
    set(fc_axes, 'ColorMap', colorMap);
    set([lr_ax, rl_ax],'CLim',[-cl cl]);
    cb = colorbar(fc_axes(end),'EastOutside');
    %cb.Layout.Tile = 'west';
    %title(fc_axes(2), app.PlotFCs.Value,'FontSize', 12);
    linkaxes(fc_axes,'xy');
   
    
    
end


%% Misc funcs

function H = matrix_polynomial(coeffs, S)
    [m, n] = size(S);
    H = zeros(n, n);
    
    for l = 1:length(coeffs)
        H = H + coeffs(l) * S^(l-1);
    end

end
