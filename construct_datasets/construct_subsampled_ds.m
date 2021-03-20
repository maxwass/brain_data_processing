%% construct datasets with different frequency filtering

include_subcortical = true;
width = 15;

if include_subcortical
    num_freqs = 87;
else
    num_freqs = 68;
end

all_intervals_to_remove ={};
all_intervals_to_remove{end+1} = {[width,num_freqs]};
low = width;
high = low+width;

while (low+width-1)<=num_freqs
    all_intervals_to_remove{end+1} = {[2,low+1],[low+width,num_freqs]};
    low = low + width;
end

all_intervals_to_remove{end+1} = {[2,num_freqs-width]};

for l = 1:length(all_intervals_to_remove)
    intervals_to_keep = all_intervals_to_remove{l};
    fc_subsample_ds(intervals_to_keep, include_subcortical)
    
end


