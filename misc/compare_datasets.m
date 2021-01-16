%% compare two datasets

close all;
clear all;

plot_write_path = '~/Documents/MATLAB/brain_data_preprocess/plots/old_vs_new_plots';
%{
%scs
% all_id (1065x1 double)
% loaded_tensor_sub (4-D double)
load('HCP_subcortical_CMData_desikan.mat')
SCs = squeeze(loaded_tensor_sub(:,:,1,:)); % (87x87x1065 double)
subject_list_sc = int64(all_id);
clear('mode3', 'loaded_bin_network_sub', 'all_id','loaded_tensor_sub'); %unused variables
%}


% OLD provided FCS
% all_fc (1x1065 cell)
% ChosenROI_cortical (1x68 double)
% ChosenROI_subcortical (1x68 double)
% subjList (1065x1 double)
load('desikan_fc_all.mat')
given_fcs  = all_fc;
subject_list_old_fc = int64(subjList);
FCs = reshape(cell2mat(all_fc),87,87,[]); %(87x87x1058 double)
clear all_fc ChosenROI_cortical ChosenROI_subcortical subjList

%make new tensor same shape as scs and fill missing entries with NaN matrix
FCs_with_nans = zeros(87,87,1065);
missing_data = [239,297,351,387,639,870,1064];
i_sub_index = 1;
for i_index = 1:1065
    fc = given_fcs(i_index);
    
    if isempty(fc{1}) %all_ismember(i_index, missing_data) %better to check if all_fc(i_index) is empty
        %put in all nana
        FCs_with_nans(:,:,i_index) = nan(87,87);
    else
        FCs_with_nans(:,:,i_index) = FCs(:,:,i_sub_index);
        i_sub_index = i_sub_index + 1;
    end
end

data_folder     = '~/Documents/MATLAB/brain_data_preprocess/data';
addpath(data_folder);
computed_fcs_given_scs =  load("data/brain_dataset_sc_fc_pairs.mat");

% NEW computed FCs
fcs_corr_mean = computed_fcs_given_scs.fcs_corr_mean;
SCs = computed_fcs_given_scs.transform_scs;
subject_list_new = computed_fcs_given_scs.final_subject_list;
subcortical_first = computed_fcs_given_scs.subcortical_first;


%new fcs are superset (i think) of old fcs. Only use ones that have old fc
for j_index = 1:length(subject_list_old_fc)

subject = subject_list_old_fc(j_index);
if ismember(j_index, missing_data) %missing from old fc
    continue
end

txt = sprintf('patient %d', subject);
fprintf('%dth patient %d...\n',j_index, subject);


%create plot with 2x3 tiles
%first row being full 87x87 sc, old fc, new fc
%second row being only cortical 68x68 sc, old fc, new fc
fig = figure('visible','off');

t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
title(t,txt, 'FontSize', 25)

nexttile
%index_in_subject_list_sc = int64(find(subject_list_sc == subject));
index_in_subject_list_new = int64(find(subject_list_new == subject));
imshow(SCs(:,:,index_in_subject_list_new)/9.9);
title('SC/9.9', 'FontSize', 15);
ylabel('Subcortical AND Cortical', 'FontSize', 15); 

nexttile
index_in_subject_list_old_fc = int64(find(subject_list_old_fc == subject));
old_corr = FCs_with_nans(:,:,index_in_subject_list_old_fc);
imshow(abs(old_corr));
title('|FC OLD Corr|','FontSize', 15);

nexttile
index_in_subject_list_new = int64(find(subject_list_new == subject));
new_corr = fcs_corr_mean(:,:,index_in_subject_list_new);
imshow(abs(new_corr));
title('|FC NEW Correlation|','FontSize', 15);

%only cortical
nexttile
%index_in_subject_list_sc = int64(find(subject_list_sc == subject));
index_in_subject_list_new = int64(find(subject_list_new == subject));
sc = SCs(:,:,index_in_subject_list_new);
sc_cortical = sc(20:end, 20:end);
imshow(sc_cortical/9.9);
%title('SC cortical (divided by 9.9)','FontSize', 15);
ylabel('ONLY Cortical', 'FontSize', 15) ;

nexttile
index_in_subject_list_old_fc = int64(find(subject_list_old_fc == subject));
old_corr_cortical = FCs_with_nans(:,:,index_in_subject_list_old_fc);
old_corr_cortical = old_corr_cortical(1:68,1:68); %see Yang script
imshow(abs(old_corr_cortical));
%title('FC cortical OLD Correlation (we only used cortical)','FontSize', 15);

nexttile
index_in_subject_list_new = int64(find(subject_list_new == subject));
new_corr_cortical = fcs_corr_mean(:,:,index_in_subject_list_new);
new_corr_cortical = new_corr_cortical(20:end, 20:end); %see processing script and 'subcortical_first'
imshow(abs(new_corr_cortical));
%title('FC cortical NEW Correlation','FontSize', 15);


%tensor of only old correlation repeated
%old_corr_cortical_repeat = repmat(old_corr_cortical, 1,1,1058);
%new_corr_cortical_repeat = 
fn = plot_write_path + "/" + string(subject) + ".jpg";
exportgraphics(t,fn,'BackgroundColor','white')

end


