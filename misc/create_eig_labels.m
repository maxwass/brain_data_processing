function [eig_labels] = create_eig_labels(eigvals, sig_digs)

%% eigenvalues can be very close together, making reading tick labels hard.
% 1) Find dense area of eigenvalues, and only label beginning and end of
% interval
% 2) shorten label length
eig_diffs = diff(eigvals);
diff_threshold = prctile(eig_diffs, 50);


diff_ = eig_diffs < diff_threshold;



start_chunk_label = find(diff_,1,'first');
end_chunk_label   = find(diff_,1,'last');
middle_chunk_label = ceil((end_chunk_label-start_chunk_label)/2);

eig_labels = {};

for i = 1:length(eigvals)
    e = eigvals(i);
    s = num2str(e); % '%0.3940'
    last_idx = min(length(s),sig_digs);
    s = s(1:last_idx);
    
    % if label is in dense region and not middle_chunk_label, clear it
    if (i > start_chunk_label) && (i < end_chunk_label) && (i~=middle_chunk_label)
        s = '';
    end

    eig_labels{i}=s;
end


end

