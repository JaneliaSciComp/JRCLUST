%--------------------------------------------------------------------------
% 2017/12/15 JJJ: Deal with overlapping spikes
function spikeFeatures = sort_overlap_(S0, S_clu, P)
    global spikeFeatures

    if ~all(size(spikeFeatures) ~= S0.featureDims)
        spikeFeatures = getSpikeFeatures(P);
    end

    % find which spikes it's overlapping with
    % only use larger unit to correct smaller unit
    % S_overlap = find_overlap_(S0.spikeTimes, S0.spikeSites, S_clu, P);
    % load raw traces and redetect from filtered dtraces. load pairs
    % Find pair of spikes by unit. For unit A, find other units that collide with spikes. which spikes are colliding with which units
    % subtract other spike waveforms by using time delays, adjust feature amplitudes, larger corecting smaller spikes, mutual corection
    t1 = tic;
    [cviSpk_o_1, cviSpk_o_12, cviDelay1] = find_overlap_(S0, S_clu, P);
    fprintf('overlapping spike detection took %0.1fs\n', toc(t1));
    %-----
    % Cancel overlap
    % subtract other spike waveforms and features
    spikeFeatures = cancel_overlap_(cviSpk_o_1, cviSpk_o_12, cviDelay1, S0, S_clu, P);
    % write to spikeFeatures
    % strategy 1. find nearest features and copy cluster identity
    % strategy 2. recluster the whole thing after cancellation

    %-----
    % recluster
    S_clu = fet2clu_(S0, P);
end % function
