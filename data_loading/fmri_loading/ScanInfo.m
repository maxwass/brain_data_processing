classdef ScanInfo
    properties
        subject_id int64 {mustBePositive}
        atlas string {mustBeMember(atlas,{'desikan'})}
        task string %{mustBeMember(task,{'REST1'})}
        scan_direction string {mustBeMember(scan_direction,{'LR', 'RL'})}
        include_subcortical logical
        rawdatafolder string = "/Volumes/Elements/brain_data"
        
    end
    methods
        function obj = ScanInfo(subject_id, atlas, task, scan_direction, include_subcortical)
            if nargin == 5
                %filter subject
                if ~isnumeric(subject_id)
                    obj.subject_id = str2double(subject_id);
                else
                    obj.subject_id = subject_id;
                end
                obj.atlas = atlas;
                %filter task
                if contains(task, 'REST1')
                    task = 'REST1';
                elseif contains(task, 'REST2')
                    task = 'REST2';
                else
                    error('%s not supported', task);
                end
                obj.task = task;
                obj.scan_direction = scan_direction;
                obj.include_subcortical = include_subcortical;
            end
        end
        function dtseries = load_functional_dtseries(obj)
            dtseries = load_functional_dtseries(obj.subject_id, obj.atlas, obj.task, obj.scan_direction, obj.rawdatafolder);
            dtseries   = dtseries(obj.get_roi_idxs(), :);
        end
        function A = extract_sc(obj)
            [A] = extract_sc(obj.subject_id, obj.atlas, obj.include_subcortical);
        end
        function [GFT, evals_vec, S] = extract_GFT(obj, GSO)
            [GFT, evals_vec, S] = extract_GFT(obj.subject_id, obj.atlas, obj.include_subcortical, GSO);
        end
        function [roi_idxs] = get_roi_idxs(obj)
            roi_idxs = get_roi_idxs(obj.atlas, obj.include_subcortical);
        end
        function [is_eq] = eq(obj, other)
            is_eq = ...
                isequal(obj.subject_id, other.subject_id) && ...
                isequal(obj.atlas, other.atlas) && ...
                isequal(obj.task, other.task) && ...
                isequal(obj.scan_direction, other.scan_direction) && ...
                isequal(obj.include_subcortical, other.include_subcortical) && ...
                isequal(obj.rawdatafolder, other.rawdatafolder);
            
            
        end
    end
end