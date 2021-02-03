function [A] = extract_sc(subject, atlas, include_subcortical)
%% load sc for subject from stored tensor file sent by Zhengwu

if isa(subject, 'char') || isa(subject,'string')
    subject = int64(str2double(subject));
end
    

%% load scs
if atlas=="desikan"
    scs_file = load('scs_desikan.mat');
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

