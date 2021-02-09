function [dist] = psd_dist(A, B, m)
% Compute distance between positive (semi?)definite matrix A and B 
%  according to distance metric m


%papers of interest
%https://proceedings.neurips.cc/paper/2014/file/f7664060cc52bc6f3d620bcedc94a4b6-Paper.pdf
%https://arxiv.org/pdf/1407.1120.pdf
%https://arxiv.org/pdf/0910.1656.pdf
    
switch lower(m)
    case 'herdin' % https://ieeexplore.ieee.org/document/1543265
        l = trace(A*B);
        p = norm(A,'fro') * norm(B,'fro');
        dist = 1 - (l/p);
    case 'airm'
        A_inv = sqrtm(inv(A));
        inter = logm(A_inv*B*A_inv);
        dist = norm(inter,'fro');
        
        %dist = norm( logm(A*inv(B)), 'fro')^2; % page 5: https://arxiv.org/pdf/1407.1120.pdf
    case 'stein' %Stein Metric
        dist = log( det( (A+B)/2 ) ) - 0.5*log( det(A*B) ); % page 5: https://arxiv.org/pdf/1407.1120.pdf
        
end

%TODO

%Log Euclidean
%Cholesky
%Power Euclidean distance


end

