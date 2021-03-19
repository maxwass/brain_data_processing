function [A] = extract_sc(subject, atlas, include_subcortical) %roi interest
%% load sc for subject from stored tensor file sent by Zhengwu

if isa(subject, 'char') || isa(subject,'string')
    subject = int64(str2double(subject));
end
    

%% load scs
if atlas=="desikan"
    
    %have we loaded the sc file (in the base workspace) yet?
    sc_file_loaded = evalin("base","exist('scs_file', 'var')");
    if sc_file_loaded
        scs_file = evalin("base","scs_file");
    end
    
    % if file dne or it exists but wrong atlas...
    if ~sc_file_loaded || ~strcmp(scs_file.atlas, atlas)
            scs_file = load('scs_desikan.mat');
        
            %put into workspace so don't have to keep reloading.
            assignin('base', 'scs_file', scs_file)
    end

elseif atlas=="destrieux"
    error("NOT IMPLIMENTED")
else
    error("NOT IMPLIMENTED")
end
    

%% find index of subject in tensor
idx = find(scs_file.subject_list==subject);
A_full_raw = scs_file.scs(:,:, idx);
A_full_transform = log(A_full_raw+A_full_raw' +1); %transform reccomended by Zhengwu

A_cortical_transform = A_full_transform(20:end, 20:end);

if include_subcortical == 1
    A = A_full_transform;
else
    A = A_cortical_transform;
end



end

