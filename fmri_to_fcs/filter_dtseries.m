function [filtered_dtseries] = filter_dtseries(dtseries, subject, atlas, include_subcortical, GSO, which_filters)
%Filter dtseries

    
if ismember("lpf", which_filters)
    %take GFT of signals, remove high frequencies, take iGFT
    [dtseries_freq, GFT, evals] = apply_GFT(dtseries, subject, atlas, include_subcortical, GSO);
    
    %extract first _ rows OR zero out last _ rows
    
    
    %apply iGFT
    filtered_dtseries = ;
    
end

if ismember("energy", which_filters)
    %remove observations with total energy above ___
    
end

if ismemmber("active_frontal_lobe", which_filters)
    % remove observations where clear thinking is occuring
    
end

filtered_dtseries = inputArg1;
end

