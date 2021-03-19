function [labels, idxs] = create_labels(vals, sig_digs, k)

%% vals can be very close together, making reading tick labels hard.
% 1) Find dense area of vals, and only label beginning and end of
% interval
% 2) shorten label length
diffs = diff(vals);
diff_threshold = prctile(diffs, 50);


diff_ = diffs < diff_threshold;



start_chunk_label = find(diff_,1,'first');
end_chunk_label   = find(diff_,1,'last');
middle_chunk_label = ceil((end_chunk_label-start_chunk_label)/2);
%k = 7;
sparse_labels = middle_chunk_label - [-2*k,-k, k, 2*k] ;

labels = {};
idxs = {};

for i = 1:length(vals)
    e = vals(i);
    e = round(e,2); %
    s = num2str(e); % '%0.3940'
    last_idx = min(length(s),sig_digs);
    s = s(1:last_idx);
    
    idx_s = num2str(i);
    
    % if label is in dense region and not middle_chunk_label, clear it
    if ismember(i, sparse_labels)
        %expection for these cases
    elseif (i > start_chunk_label) && (i < end_chunk_label) && (i~=middle_chunk_label)
        s = '';
        idx_s = '';
    end
    
    %make expection for middle_chunk_label +/- 5
    

    labels{i}=s;
    idxs{i} = num2str(idx_s);
end


end

