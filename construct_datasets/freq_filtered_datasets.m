%% construct datasets with different frequency filtering

include_subcortical = true;
width = 15;

if include_subcortical
    num_freqs = 87;
else
    num_freqs = 68;
end

all_intervals ={};
all_intervals{end+1} = {[width,num_freqs]};
low = width;
high = low+width;

while (low+width-1)<=num_freqs
    all_intervals{end+1} = {[2,low+1],[low+width,num_freqs]};
    low = low + width;
end

all_intervals{end+1} = {[2,num_freqs-width]};

for l = 1:length(all_intervals)
    intervals_to_keep = all_intervals{l};
    fc_subsample_ds(intervals_to_keep, include_subcortical)
    
end


