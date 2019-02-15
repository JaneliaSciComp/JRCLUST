function [sites1, sites2, sites3] = getSecondaryPeaks(obj)
    %GETSECONDARYPEAKS
    minVals = squeeze(min(obj.meanWfGlobal) - obj.meanWfGlobal(1, :, :));

    [~, sites1] = min(minVals);

    minVals(sub2ind(size(minVals), sites1, 1:numel(sites1))) = 0;
    [~, sites2] = min(minVals);

    minVals(sub2ind(size(minVals), sites2, 1:numel(sites2))) = 0;
    [~, sites3] = min(minVals);
end