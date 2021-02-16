function [intervals] = parse_freq_range(raw_text)
%text of format [low-high],[low2-high2], ...,[lowk,highk]
%cell array pf  {[low,high], ...

%cast to char vector
text = char(raw_text);

%remove whitespace
text(text==' ') = [];

%basic testing
if ( text(1)~='[' ) || ( text(end)~=']' )
    error('freq intervals incorrect: %s\n. Format: [1-6],[7,10]', raw_text)
end

%parse
text_intervals = split(text,',');
intervals = cell(1,length(text_intervals));

%check that intervals are not intersecting. Assumes intervals come in sorted
%order
check = zeros(1,2*length(text_intervals));


old_low = -1;
old_high = -1;
for i = 1:length(text_intervals)
    a = text_intervals{i};
    a(1)  = [];
    a(end)=[];
    low_and_high = split(a,'-');
    if length(low_and_high)~=2
        error('splitting %s into more than 2 numbers!', a);
    end
    low  = str2double(low_and_high{1});
    high = str2double(low_and_high{2});
    if (low>=high)
        error('invalid interval: low(%d) > high(%d)',low, high);
    elseif (low <= old_high)
        error('invalid intervals: cannot intersect. low(%d) <= old_high(%d)')
    end
    intervals{i} = [low, high];
    old_high = high;
end
    
   
    
end

