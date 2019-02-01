function weightedLoc = getWeightedLoc(featureWeights, positions)
    %GETWEIGHTEDLOC Get feature-weighted location of clusters from positions
    if isrow(positions)
        positions = positions';
    end

    featureWeights = featureWeights.^2;
    sos = sum(featureWeights);

    weightedLoc = sum(bsxfun(@times, featureWeights, positions))./sos;
end