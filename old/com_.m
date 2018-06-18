%--------------------------------------------------------------------------
function vrCentroid = com_(vrVpp, vrPos)
    vrVpp_sq = vrVpp(:).^2;
    vrCentroid = sum(vrVpp_sq .* vrPos(:)) ./ sum(vrVpp_sq);
end %func
