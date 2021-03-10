function [eigval_idxs_to_remove] = str_to_logical_vec_freqs(freq_remove_str, eigvals)
% freq_remove_str :: string of the form '[-100.2, -9.8], [-8,100]'
% eigvals :: vector of increasing eigenvalues

% return logical vector indicating which eigenvalue indices to keep

    % remove final ']'
	freq_remove_str = freq_remove_str(1:end-1);
    
    % split by ']'
	interval_string_cells = strsplit(freq_remove_str,']');
    
    % find index of '[' and remove everything before '['
    interval_string_cells = cellfun(@(s) extractAfter(s,'['), interval_string_cells, 'UniformOutput', false);
    interval_string_cells = cellfun(@(s) erase(s,' '), interval_string_cells, 'UniformOutput', false);
    
    % split numbers by ','
	intervals_string_cells_sep = cellfun(@(s) strsplit(s,','), interval_string_cells, 'UniformOutput', false);
	
    % get rid of excess ','
    %interval_string_cells_sep = cellfun(@(s) erase(s,','), interval_string_cells, 'UniformOutput', false);
    
   
    % convert into array of doubles at each cell location
	freq_intervals_to_remove = cellfun(@(c) [str2double(c{1}), str2double(c{2})], intervals_string_cells_sep, 'UniformOutput', false);
	
    % sanity check.
    %  -intervals must be disjoint. NOT CHECKED YET 
    %  -must be valid intervals: low < high
    for l = 1:length(freq_intervals_to_remove)
        interval = freq_intervals_to_remove{l};
        if interval(1) >= interval(2)
            error('Given interval in which low >= high. Must increase.')
        end
    end
    
    
    fprintf('removing freqs in intervals...\n')
	for l = 1:length(freq_intervals_to_remove)
        interval = freq_intervals_to_remove{l};
        fprintf('[%f, %f] ', interval(1), interval(2));
    end
	fprintf('\n')
                
	%convert  cell array of intervals => logical vector over idxs
	eigval_idxs_to_remove = false(size(eigvals));
	for l = 1:length(freq_intervals_to_remove)
        low = freq_intervals_to_remove{l}(1);
        high = freq_intervals_to_remove{l}(2);
        new_to_remove = (low <= eigvals) & (eigvals <= high);
        eigval_idxs_to_remove = new_to_remove | eigval_idxs_to_remove;
    end
    
    eigval_idxs_to_keep = ~eigval_idxs_to_remove;
    
end

