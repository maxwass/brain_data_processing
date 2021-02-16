%% Explore and Viz Brain Data
clear all; close all;
N_full = 87; N_sub = 68;
% loads 6 matrices: 
% func_corr_full = 1058 functional correlation matrices. Each row 
%                   is reshaped 87x87 matrix.
% func_corr_sub = 1058 functional correlation matrices. Each row 
%                   is reshaped 68x68 matrix. Subset of above.
% struct_full = 1058 structual connectivity matrices. Each row 
%                   is reshaped 87x87 matrix.
% struct_sub = 1058  structual connectivity matrices. Each row 
%                   is reshaped 68x68 matrix. Subset of above.
%struct_sparsity_full = the sparsity of all 1085 patient full strctural
%                       connectivty graphs
%struct_sparsity_sub = the sparsity of all 1085 patient sub strctural
%                       connectivty graphs

% to use, pull out a row and reshape it to square size.

load brain_data.mat

i = 250;
sc_sub_row = struct_sub(i,:);
sc_sub     = reshape(sc_sub_row, N_sub, N_sub);


%% sparsity of structural graphs
[max_val, max_ind] = max(struct_sparsity_sub); %index 574
[min_val, min_ind] = min(struct_sparsity_sub); %index 324
edges   = linspace(min_val, max_val, 200);
figure();
histogram(struct_sparsity_sub,edges);
ax = gca;
ax.YAxis.FontSize = 15; %for y-axis
ax.XAxis.FontSize = 15; %for y-axis
xlabel("Sparsity: # edges/2 / N(N-1)/2", 'FontSize', 25);
ylabel("Number of Patients", 'FontSize', 25);
title("Distribution of Structural Connectivity Graph Sparsities", 'FontSize', 30);



%% Distribution of egde weights for patient i
total_poss_edges      = N_sub*(N_sub-1)/2;
nonzero_indx          = (triu(sc_sub,1)>0);
non_zero_edge_weights = sc_sub(nonzero_indx);
num_zeros_edges       = total_poss_edges - length(non_zero_edge_weights);
all_edges             = cat(1,zeros(num_zeros_edges,1),non_zero_edge_weights);
    
[max_edge_weight, ~] = max(all_edges);
bins   = linspace(0, max_edge_weight, 20);
figure();
histogram(all_edges,bins);
%histogram(non_zero_edge_weights,bins);
t = sprintf(['Distrib of Edge Weights in SC Graph: Patient %i, sparsity: ',num2str(struct_sparsity_sub(i))], i);
title(t, 'FontSize', 30);
xlabel("Edge Weights", 'FontSize', 25);
ylabel("Number of Edges", 'FontSize', 25);
pause(0.1);
   
    
%% slowly remove small edge weights
    
figure();
t = sprintf('Slowly Remove Small Edges of Structural Connectivity Graph: Patient %i', i);
sgtitle(t, 'FontSize', 30);
cutoffs = quantile(non_zero_edge_weights,4); 
    
subplot(5,1,1);
adj = sc_sub;
G = graph(adj,'omitselfloops');
p = plot(G,'Layout','force');%'force');%,'WeightEffect','direct');
title('all edges', 'FontSize', 20)
%layout(p,'subspace','Dimension',50);

subplot(5,1,2);
adj(adj<cutoffs(1)) = 0;
G = graph(adj,'omitselfloops');
p = plot(G,'Layout','force');
%p = plot(G,'Layout','circle');
title(sprintf('edges above 20th quantile <==> edges above %f',cutoffs(1)),'FontSize', 20);
    
subplot(5,1,3);
adj(adj<cutoffs(2)) = 0;
G = graph(adj,'omitselfloops');
p = plot(G,'Layout','force');
%p = plot(G,'Layout','circle');
title(sprintf('edges below 40th quantile <==> edges above %f',cutoffs(2)),'FontSize', 20);
   
subplot(5,1,4);
adj(adj<cutoffs(3)) = 0;
G = graph(adj,'omitselfloops');
p = plot(G,'Layout','force');
%p = plot(G,'Layout','circle');
title(sprintf('edges below 60th quantile <==> edges above %f',cutoffs(3)),'FontSize', 20);
    
subplot(5,1,5);
adj(adj<cutoffs(4)) = 0;
G = graph(adj,'omitselfloops');
p = plot(G,'Layout','force');
%p = plot(G,'Layout','circle');
title(sprintf('edges below 80th quantile <==> edges above %f',cutoffs(4)),'FontSize', 20);
    
%p = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight);
%layout(p,'force','WeightEffect','direct')
%layout(p,'force3','Iterations',10)

