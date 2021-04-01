%% Sanity checking SC graphs. 
% Do any contain NaNs? Any all zeros? Any with negative edges? Any nodes
% have no edges? Any disconnected?


% subject_list (1065x1 int64)
% scs (87x87x1065 double)
sc_file = load('data/scs_desikan.mat');
SCs_upper_tri = sc_file.scs; 
subject_list_sc = sc_file.subject_list;

SCs = apply_to_tensor_slices(@(x) (x+x'), SCs_upper_tri);
cortical_idxs = get_roi_idxs(atlas, include_subcortical);
cortical_SCs = SCs(cortical_idxs,cortical_idxs,:);
SCs_transform = apply_to_tensor_slices(@(x) log(x+x'+1), SCs_upper_tri);


%% any graphs have nan?
anynan = any(any(any(isnan(SCs))));
assert(~anynan, 'some SC contains NaN values');

%% any matrix have all zeros?
allnonzero = all(any(any(SCs>0,2)));
assert(allnonzero, 'some SC consists of all zeros');

%% any matrix have negative edges?
anyneg = any(any(any(SCs<0,2)));
assert(~anyneg, 'some SC contains negative values');


%% Any total graphs disconnected? => 35
which_total_connected = apply_to_tensor_slices(@check_connectivity, SCs);
num_disconnected_total = sum(~which_total_connected);

%% Any cortical graphs disconnected? => 0
which_cortical_connected = apply_to_tensor_slices(@check_connectivity, cortical_SCs);
num_disconnected_cortical = sum(~which_cortical_connected);
assert(num_disconnected_cortical==0, 'some CORTICAL SC(s) are disconnected');


%% any nodes in total graphs with no neighbors? => 35 patients, all Node 8
[patients_with_zero_degree, nodes_with_zero_edges] = zero_degree_patients(SCs, subject_list_sc);
display_patient_with_zd_nodes(patients_with_zero_degree, nodes_with_zero_edges, 'TOTAL');

%% any nodes in CORTICAL graphs with no neighbors? => 0 patients
[patients_with_zero_degree, nodes_with_zero_edges] = zero_degree_patients(cortical_SCs, subject_list_sc);
display_patient_with_zd_nodes(patients_with_zero_degree, nodes_with_zero_edges, 'CORTICAL');
assert(isempty(patients_with_zero_degree), 'some CORTICAL SC(s) have nodes with 0 edges');


function display_patient_with_zd_nodes(patients_with_zero_degree, nodes_with_zero_edges, txt)
    if isempty(patients_with_zero_degree)
        return
    end
    fprintf('These %d patients have at least one node with zero degree in %s graph\n', length(patients_with_zero_degree), txt);
    for p = 1:length(patients_with_zero_degree)
        fprintf('patient %d: nodes with zero degree', patients_with_zero_degree(p));
        disp(cell2mat(nodes_with_zero_edges(p)));
    end
    fprintf('\n');
end
function [patients_with_zero_degree, nodes_with_zero_edges] = zero_degree_patients(SCs_tensor, subject_list_sc)
    degrees = sum(SCs_tensor, 2);
    zero_degree_slice = @(d, eps) any( abs(d-0.0)<eps );
    scs_with_zero_degree = zero_degree_slice(degrees, .00001);
    zero_degree_idxs = find(scs_with_zero_degree);
    patients_with_zero_degree = subject_list_sc(zero_degree_idxs);
    
    nodes_with_zero_edges = cell(length(patients_with_zero_degree),1);
    %find the nodes which have 0 degree
    for l = 1:length(zero_degree_idxs)
        idx = zero_degree_idxs(l);
        ds = degrees(:,:,idx);
        nodes_with_zero_edges{l} = find( abs(ds)<.00001);
    end
end
