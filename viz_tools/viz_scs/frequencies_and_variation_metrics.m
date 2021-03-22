%% for each SC, compute the eigenvalues and corresponding 
%  Zero-Crossings/Total Variations for each frequency.
clear; close;

atlas = 'desikan';
GSOs = ["A","A_norm", "L", "L_norm"];
%include_subcortical_list = [false, true];
include_subcortical = false;
scs_file = load('scs_desikan.mat');
subject_list = scs_file.subject_list;


for idx = 1:length(subject_list)
    subject_id = subject_list(idx);
    gso_struct = struct(GSOs(1),[], GSOs(2), [], GSOs(3), [], GSOs(4), []);
    for gso_idx = 1:length(GSOs)
        s = struct('eigenvals', [], 'total_variations', [], 'zero_crossings', []);
        
        GSO = GSOs(gso_idx);
        [GFT, s.eigenvals, S] = extract_GFT(subject_id, atlas, include_subcortical, GSO);
        A = extract_sc(subject_id, atlas, include_subcortical);

        %% interpretation of freq components: zero crossings and total variation
        V = GFT';
        L_sign = diag(sum(sign(A),2))-sign(A);
        % ith position is # 0-crossings of ith freq component
        s.zero_crossings = zero_crossings(V, A);

        L = diag(sum(A,2)) - A;
        s.total_variations = total_variation(V, L);
        
        % place in larger struct
        gso_struct.(GSO) = s;
        %preview_freq_plot(GSO, s.eigenvals, s.zero_crossings, s.total_variations);
        %close;
    end
    
    % save
    filepath = fullfile("viz_scs", "scs_spectral_info", string(subject_id));
    save(filepath, 'gso_struct', 'atlas', 'include_subcortical');
    fprintf('%d/1065: %d\n', idx, subject_id);
    
end


%% Misc funcs

function preview_freq_plot(GSO, eigvals, zero_crossings, total_variation)
    figure;
    plot(eigvals, zero_crossings, 'DisplayName', 'Zero-Crossings');
    hold on;
    plot(eigvals, total_variation, 'DisplayName', 'Total Variation');
    hold off;
    xlabel('Eigenvalues');
    ylabel('Measure of Variation');
    title(GSO);

end
