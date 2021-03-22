classdef ScanInfo
    properties
        subject_id int64 {mustBePositive}
        atlas string {mustBeMember(atlas,{'desikan'})}
        task string {mustBeMember(task,{'REST1'})}
        scan_direction string {mustBeMember(scan_direction,{'LR', 'RL'})}
        include_subcortical logical
        rawdatafolder string = "/Volumes/Elements/brain_data"
        
    end
    methods
        function obj = ScanInfo(subject_id, atlas, task, scan_direction, include_subcortical)
            if nargin == 5
                obj.subject_id = subject_id;
                obj.atlas = atlas;
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
    end
end