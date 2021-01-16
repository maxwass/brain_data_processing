function [sc_full, sc_cortical] = transform_sc(sc)
% Take raw structural connectivity readings from preprocessing
%pipeline (and upper triangular matrix with values in (0, 10000??) and
%output a symmetric matrix with transformed values
%   Detailed explanation goes here

%sc is an upper triangular matrix
% apply transformation to make symmetric and decrease variation between
% small and large entries.
sc_full        = log(sc+sc'+1);
sc_cortical    = sc_full(20:end,20:end);

end

