%% fix cached data to reflect new chosen_roi fileds

clear; clc;

cached_dir = 'data/cached_desikan/rfMRI_REST1';


list_file_structs = dir('data/cached_desikan/rfMRI_REST1/*_*.mat');

for i = 1:length(list_file_structs)
    start = tic;
    fname = list_file_structs(i).name;
    fn = [list_file_structs(i).folder '/' fname];
    fdata = load(fn);
    
    if ~isfield(fdata, 'subcortical_first')
        fprintf('%d: %s already processed\n', i, fname);
        continue;
    end
        
    %fix chosen roi
    chosen_roi = struct('cortical', fdata.chosen_roi.cortical, ...
        'subcortical', fdata.chosen_roi.subcortical,...
        'subcortical_first', logical(fdata.subcortical_first));
    
    
    atlas = fdata.atlas;
    scan = fdata.scan;
    tasktype = fdata.tasktype;
    dtseries = fdata.dtseries;
    
    fn = [list_file_structs(i).folder '/' fname];
    save(fn, 'atlas','scan', 'tasktype','chosen_roi', 'dtseries');
    
    stop = toc(start);
    fprintf('%d: %s fixed in %.2f (s)!\n', i, fname, stop);
end
