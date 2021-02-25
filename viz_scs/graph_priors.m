%% construct 'prior' of graphs
% which edges do all/most patients have, and what is the distribution and
% summary statistics of those values?
clear;clc;close all;

scs = load('scs_desikan.mat').scs;

%% find the edges that many scs share
fraction = .8; %what fraction of scs should have this?
common_edges_upper_tri = sum((scs>0),3) >= fraction*length(scs);
common_edges  = common_edges_upper_tri + common_edges_upper_tri';

%% find summary statistic of all edge distributions
edge_medians = transform(median(scs,3));
edge_means   = transform(mean(scs,3));
edge_max     = transform(max(scs,[],3));
edge_min     = transform(min(scs,[],3));

%% choose one of the summary statistics
edge_values = edge_medians;

%% construct 'prior'
prior = common_edges.*edge_values;



%% plotting

f = figure();
t = tiledlayout('flow');
prior_axes = gobjects(5,1); 

%inspect distribution of edge weights on between node p and q
p = 30;
q = 20;
p_q_edges = get_edge_weights_list(p,q,scs);
ax = nexttile();
histogram(p_q_edges,'Normalization','probability');
title(sprintf('Edge Weights between %d and %d', p,q));


l=1;
prior_axes(l) = nexttile();
imagesc(common_edges);
colormap( prior_axes(l), [0,0,0; 1,1,1] )
txt = sprintf("White implies >=%.1f%s of patients have nonzero edge", fraction*100, '%');
title(txt);


l=l+1;
prior_axes(l) = nexttile();
imagesc(common_edges.*edge_medians);
title('Medians');

l=l+1;
prior_axes(l) = nexttile();
imagesc(common_edges.*edge_means);
title('Means');

l=l+1;
prior_axes(l) = nexttile();
imagesc(common_edges.*edge_max);
title('Max');

l=l+1;
prior_axes(l) = nexttile();
imagesc(common_edges.*edge_min);
title('Min');

cl = max(max(edge_max));
set(prior_axes(2:end),'CLim',[0 cl]);
cb = colorbar(prior_axes(end),'EastOutside');
linkaxes(prior_axes,'xy');
for k = 1:length(prior_axes) %make square
    daspect(prior_axes(k),[1 1 1]);
end

% Given summary statistic is upper triangular and has large values. 
%   Make symmtric and squash values
function B = transform(A)
    B = log(A' + A);
end

function edges = get_edge_weights_list(p,q,scs)
    if q<=p
        edges = log(scs(q,p,:)+1);
    else
        edges = log(scs(p,q,:)+1);
    end
end