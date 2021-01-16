%% compare two datasets
path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo)); %recursively adds all repo files/folders

close all;
clear all;

plot_write_path = '~/Documents/MATLAB/brain_data_preprocess/plots/old_vs_new_plots';
%{
%scs
% all_id (1065x1 double)
% loaded_tensor_sub (4-D double)
load('scs.mat')
SCs = squeeze(loaded_tensor_sub(:,:,1,:)); % (87x87x1065 double)
subject_list_sc = int64(all_id);
clear('mode3', 'loaded_bin_network_sub', 'all_id','loaded_tensor_sub'); %unused variables
%}


% OLD provided FCS
old_fc_file = load('data/correlations_desikan_old.mat');
%given_fcs  = all_fc;
subject_list_old_fc = int64(old_fc_file.subject_list);
FCs = old_fc_file.fcs; %(87x87x1058 double)
clear old_fc_file


computed_fcs_given_scs =  load("data/brain_dataset_sc_fc_pairs.mat");

% NEW computed FCs
fcs_corr_mean = computed_fcs_given_scs.fcs_corr_mean;
SCs = computed_fcs_given_scs.transform_scs;
subject_list_new = computed_fcs_given_scs.final_subject_list;
subcortical_first = computed_fcs_given_scs.subcortical_first;


%new fcs are superset (i think) of old fcs. Only use ones that have old fc
for j_index = 1:length(subject_list_old_fc)

subject = subject_list_old_fc(j_index);

txt = sprintf('patient %d', subject);
fprintf('%dth patient %d...\n',j_index, subject);


%create plot with 2x3 tiles
%first row being full 87x87 sc, old fc, new fc
%second row being only cortical 68x68 sc, old fc, new fc
fig = figure();%('visible','off');

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
old_corr = FCs(:,:,index_in_subject_list_old_fc);
imshow(abs(old_corr));
title('|FC OLD Corr| -lr,rl,mean?','FontSize', 15);

nexttile
index_in_subject_list_new = int64(find(subject_list_new == subject));
new_corr = fcs_corr_mean(:,:,index_in_subject_list_new);
imshow(abs(new_corr));
title('|FC NEW Correlation| - mean','FontSize', 15);

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
old_corr_cortical = FCs(:,:,index_in_subject_list_old_fc);
old_corr_cortical = old_corr_cortical(1:68,1:68); %see Yang script
imshow(abs(old_corr_cortical));
%title('FC cortical OLD Correlation (we only used cortical)','FontSize', 15);

nexttile
index_in_subject_list_new = int64(find(subject_list_new == subject));
new_corr_cortical = fcs_corr_mean(:,:,index_in_subject_list_new);
new_corr_cortical = new_corr_cortical(20:end, 20:end); %see processing script and 'subcortical_first'
imshow(abs(new_corr_cortical));
%title('FC cortical NEW Correlation','FontSize', 15);


% add color bar. All values between 0/1



fn = plot_write_path + "/" + string(subject) + ".jpg";
exportgraphics(t,fn,'BackgroundColor','white')

end


