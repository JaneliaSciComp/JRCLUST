%--------------------------------------------------------------------------
function nearestSpikes = spikesNearMidpoint(spikes, spikeTimes, proportion)
    % get proportion of given spikes closest to the
    % midpoint; i.e., if nTime_clu == 4 and nSpikes == 4000, get the
    % 4000/4 == 1000 centermost spikes

    if proportion <= 0 || proportion > 1
        proportion = 1;
    end

    midpoint = round(numel(spikeTimes)/2);

    % rank spikes by distance from midpoint, favoring nearer (lower is better)
    ranking = rankOrder(abs(spikes - midpoint), 'ascend');
    maxRanking = round(numel(spikes)*proportion);
    nearestSpikes = spikes(ranking <= maxRanking);
end % function
