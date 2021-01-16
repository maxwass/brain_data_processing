%% find repeats in subject list

%there was one patient repeated due to a bug in the hcp_1200 subject list
%creation function. This found the repeat.

% Get a list of all files and folders in this folder.
files = dir('/Volumes/Elements/brain_data');
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.
fprintf('number of directories %d\n', length(subFolders));
dir_list = zeros(length(subFolders),1);
for k = 1 : length(subFolders)
  if (subFolders(k).name == ".") || (subFolders(k).name == "..")
      continue
  end
  dir_id = double(string(subFolders(k).name));
  fprintf('Sub folder #%d = %s =?= %d \n', k, subFolders(k).name, dir_id);
  if ismember(dir_id, dir_list)
      fprintf('FOUND DUPLICATE DIR: %d at %dth index \n', dir_id, k);
  end
  dir_list(k) = dir_id;
end