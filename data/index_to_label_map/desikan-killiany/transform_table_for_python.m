%% table to struct of arrays

%load in data with Import Data tool in Home
desikan_killiany_cortical_index_to_label_struct.freesurfer_label_id = ...
    desikankillianycorticalindextolabel.freesurfer_label_id;

% string array => cell array b/c python loads cell arrays as tuples
% (actualy numpy array)
desikan_killiany_cortical_index_to_label_struct.label_name = ...
    cellstr(desikankillianycorticalindextolabel.label_name);
desikan_killiany_cortical_index_to_label_struct.lobe = ...
    cellstr(desikankillianycorticalindextolabel.lobe);
save('data/index_to_label_map/desikan-killiany/desikan_killiany_cortical_index_to_label.mat', 'desikan_killiany_cortical_index_to_label_struct')