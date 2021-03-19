function [segmentation_atlas] = load_atlas(atlas, subject, rawdatafolder)
%% load the subject specific label (atlas)
% (note: here LR stands for left and right hemispheres, NOT phase
% encoding - aka directionality of scanning left to right or right to
% left.)

% rawdatafolder: str path to all brain data (atlas', fmri's, etc). Currently
%               on external hd
startload = tic;
if(strcmp(atlas,"desikan"))
    subjlab = ft_read_cifti([rawdatafolder '/' subject '/' subject '.aparc.32k_fs_LR.dlabel.nii']);
    segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc']);
elseif(strcmp(atlas,"destrieux"))
    subjlab = ft_read_cifti([rawdatafolder '/' subject '/' subject '.aparc.a2009s.32k_fs_LR.dlabel.nii']);
    segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc_a2009s']);
else
    error("Atlas " + atlas + " not found. Use Desikan or destrieux.")
end

endload = toc(startload);
fprintf('load_atlas %s - %.2f\n', subject, endload);
end

