function [logical_vec] = intervals_to_logical_vec(intervals, y)
%convert  cell array of intervals => logical vector over idxs
% y :: vector of values to which the intervals refer

% If logical_vec(i) = true <==> y(i) in at least one interval

	logical_vec = false(size(y));
    for l = 1:length(intervals)
        low = intervals{l}(1);
        high = intervals{l}(2);
        
        % are any y's in this interval?
        new_y_logical = (low <= y) & (y <= high);
        
        % if so, their corresponding index in idxs_to_remove is now true
        logical_vec = new_y_logical | logical_vec;
    end

end