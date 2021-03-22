function [intervals] = str_to_intervals(remove_str)
%Takes a string of the form: '[-100.2, -9.8], [-8,100]' and constructs a
%cell arry of intervals { {[-100.2, -9.8]}, {[-8,100]} }

    [intervals] = parse_input_str(remove_str);
    error_checking(intervals);

end

function [cell_arr_intervals] = parse_input_str(input_str)
    % remove final ']'
	input_str = input_str(1:end-1);
    
    % split by ']'
	interval_string_cells = strsplit(input_str,']');
    
    % find index of '[' and remove everything before '['
    interval_string_cells = cellfun(@(s) extractAfter(s,'['), interval_string_cells, 'UniformOutput', false);
    interval_string_cells = cellfun(@(s) erase(s,' '), interval_string_cells, 'UniformOutput', false);
    
    % split numbers by ','
	intervals_string_cells_sep = cellfun(@(s) strsplit(s,','), interval_string_cells, 'UniformOutput', false);
	
    % get rid of excess ','
    %interval_string_cells_sep = cellfun(@(s) erase(s,','), interval_string_cells, 'UniformOutput', false);
    
    % convert into array of doubles at each cell location
	freq_intervals_to_remove = cellfun(@(c) [str2double(c{1}), str2double(c{2})], intervals_string_cells_sep, 'UniformOutput', false);

    cell_arr_intervals = freq_intervals_to_remove;
end

function error_checking(cell_arr_intervals)
    % sanity check.
        %  -intervals must be disjoint. NOT CHECKED YET 
    %  -must be valid intervals: low < high
    for l = 1:length(cell_arr_intervals)
        interval = cell_arr_intervals{l};
        if interval(1) >= interval(2)
            error('Given interval in which low >= high. Must increase.')
        end
    end
end
