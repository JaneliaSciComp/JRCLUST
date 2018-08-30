%--------------------------------------------------------------------------
% 10/27/17: Detailed error report
function flag = clusterDataConsistent(S_clu)
    % Check the clustering by making sure the number of clusters is self-consistent. A canary.

    flag = 0; % assume invalid by default
    nClusters = S_clu.nClusters;

    fieldNames = fieldnames(S_clu);
    if isempty(fieldNames)
        return;
    end

    fieldNames1Dim = structFieldsByDim('S_clu', 1);
    fieldNames2Dim = structFieldsByDim('S_clu', 2);
    fieldNames3Dim = structFieldsByDim('S_clu', 3);

    % check that all fields of S_clu have nClusters elements in their respective dimensions
    errors1Dim = zeros(size(fieldNames1Dim));
    for i = 1:numel(fieldNames1Dim)
        fn = fieldNames1Dim{i};
        if isfield(S_clu, fn)
            errors1Dim(i) = numel(S_clu.(fn)) ~= nClusters;
        end
    end

    errors2Dim = zeros(size(fieldNames2Dim));
    for i = 1:numel(fieldNames2Dim)
        fn = fieldNames2Dim{i};
        if isfield(S_clu, fn)
            errors2Dim(i) = ~all(size(S_clu.(fn)) == [nClusters, nClusters]);
        end
    end

    errors3Dim = zeros(size(fieldNames3Dim));
    for i = 1:numel(fieldNames3Dim)
        fn = fieldNames3Dim{i};
        if isfield(S_clu, fn)
            errors3Dim(i) = size(S_clu.(fn), 3) ~= nClusters;
        end
    end

    flag = ~any(errors1Dim) && ~any(errors2Dim) && ~any(errors3Dim);
end % function
