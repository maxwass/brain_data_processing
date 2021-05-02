function [segmentation_atlas] = load_atlas(atlas, subject, rawdatafolder) %scan_direction
%% load the subject specific label (atlas)
% (note: here LR stands for left and right hemispheres, NOT phase
% encoding - aka directionality of scanning left to right or right to
% left.)

% rawdatafolder: str path to all brain data (atlas', fmri's, etc). Currently
%               on external hd


% must all be char's for [ ] to produce char vector, not list of strings
rawdatafolder  = char(rawdatafolder);
atlas          = char(atlas);
%scan_direction = char(scan_direction);
if isnumeric(subject)
    subject = num2str(subject);
end

startload = tic;


if(strcmp(atlas,"desikan"))
    %path2atlas = [rawdatafolder '/' subject '/' subject '.aparc.32k_fs_' scan_direction, '.dlabel.nii'];
    %path2atlas = [rawdatafolder '/' subject '/' subject '.aparc.32k_fs_LR.dlabel.nii'];
    filename = strcat(subject,'.aparc.32k_fs_LR.dlabel.nii');
    path2atlas = fullfile(rawdatafolder, subject, filename);
    dne_exception(path2atlas);
    subjlab = ft_read_cifti(path2atlas); %[rawdatafolder '/' subject '/' subject '.aparc.32k_fs_LR.dlabel.nii']); % BUG???? Hardcoded scan dir as LR!!
    load_atlas_cmd = strcat('subjlab.x',subject,'_aparc');
    segmentation_atlas = eval(load_atlas_cmd); %eval(['subjlab.x' subject '_aparc']);
elseif(strcmp(atlas,"destrieux"))
    error('%s not supported yet', atlas);
    %{
    path2atlas = [rawdatafolder '/' subject '/' subject '.aparc.a2009s.32k_fs_' scan_direction '.dlabel.nii'];
    dne_exception(path2atlas);
    subjlab = ft_read_cifti(path2atlas); %[rawdatafolder '/' subject '/' subject '.aparc.a2009s.32k_fs_LR.dlabel.nii']);
    segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc_a2009s']);
    %}
else
    error("Atlas " + atlas + " not found. Use Desikan or destrieux.")
end

endload = toc(startload);
%fprintf('load_atlas %s - %.2f\n', subject, endload);
end

function dne_exception(path2atlas)
    if ~isfile(path2atlas)
        throw(MException('Loading:DoesNotExist:Atlas', path2atlas));
    end
end

