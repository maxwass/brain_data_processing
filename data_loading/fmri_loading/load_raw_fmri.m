function [fmri_data_struct] = load_raw_fmri(path2fmri)
%takes path to raw fmri file, return fmri data
% example:
%  <-rawdatafolder-----------><subject><tasktype><scan_dir>
% '/Volumes/Elements/brain_data/996782/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'

dne_exception(path2fmri);

fmri_data_struct = ft_read_cifti(path2fmri);

end

function dne_exception(path2fmri)
    if ~isfile(path2fmri)
        throw(MException('Loading:DoesNotExist:fMRI', path2fmri));
    end
end
