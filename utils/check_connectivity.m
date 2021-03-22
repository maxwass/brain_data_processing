function isConnected = check_connectivity(A)
    g = digraph(A);
    bins = conncomp(g, 'Type', 'weak');
    isConnected = all(bins == 1);
end