function plot_fc_metrics(ax, fcs, optional)



    %max_eigs is the maximum magnitude eigenvector for each covarianace matrix
    %max_eig  = @(c) eigs(c,1);
    max_two_eigs  = @(c) eigs(c,2);
    max_eigs = apply_to_tensor_slices(max_two_eigs, fcs);

    % *_diffs is a vector of distance between adjacent matrices
    airm_diffs   = zeros(num_windows-1, 1);
    stein_diffs  = zeros(num_windows-1, 1);
    herdin_diffs = zeros(num_windows-1, 1);


    for j = 1:(num_windows-1)
        %airm_diffs(j)   = psd_dist( covs(:,:,j), covs(:,:,j+1), 'airm');
        stein_diffs(j)  = psd_dist( fcs(:,:,j), fcs(:,:,j+1), 'stein');
        herdin_diffs(j) = psd_dist( fcs(:,:,j), fcs(:,:,j+1), 'herdin');
    end
    yyaxis(ax,'left'); %main metric on left
    plot(ax, max_eigs(1,:), 'k', 'LineWidth', 1, 'DisplayName', '1 max |eigs|');
    hold(ax,'on');
    plot(ax, max_eigs(2,:), 'k:','LineWidth', 1, 'DisplayName', '2 max |eigs|');
    ylabel(ax, 'Max Eigs', 'FontSize', 15);
    set(ax, 'YColor', 'k')
    
    %plot horizontal line where cutoff is
    %place txt which shows percentile and corresponding raw value on this line
    txt = sprintf('%.0f', optional.raw_cutoff);
    if optional.use_percentile
        txt = sprintf('%s - percentile %.2f', txt, optional.percentile_val);
    end
    xline(ax, optional.raw_cutoff, 'y', txt, 'FontSize', 10,...
        'LabelHorizontalAlignment', 'center');

    %highlight xticks that were/would be removed!


    
    
    
    %{
    yyaxis(ax,'right');
    %plot(ax, [stein_diffs;mean(stein_diffs)], 'LineWidth', 1, 'DisplayName', 'stein');
    x_diffs = 1.5:1:(num_windows-.5);
    plot(ax, x_diffs, herdin_diffs, 'm', 'LineWidth', 1, 'DisplayName', 'herdin');
    ylabel(ax, 'Covariance Measures/Distances', 'FontSize', 15);
    set(ax, 'YColor', 'm') % set y axis color to same as line color
    %}
    hold(ax, 'off');
end

xlabel(ax, 'Windows')
%ylim(ax, [0,100])
legend(ax)




end

