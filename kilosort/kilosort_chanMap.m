%--------------------------------------------------------------------------
function S_chanMap = kilosort_chanMap(P)
    [fpath, ~, ~] = fileparts(P.vcFile);
    chanMap = P.chanMap;
    chanMap = chanMap(:)';
    % chanMap = [33 34 8 10 12 14 16 18 20 22 24 26 28 30 32 ...
    %     7 9 11 13 15 17 19 21 23 25 27 29 31 1 2 3 4 5 6];

    % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
    % Now we declare which channels are "connected" in this normal ordering,
    % meaning not dead or used for non-ephys data

    % connected = true(34, 1); connected(1:2) = 0;
    connected = true(size(chanMap));
    connected(P.viSiteZero) = 0;

    % now we define the horizontal (x) and vertical (y) coordinates of these
    % 34 channels. For dead or nonephys channels the values won't matter. Again
    % I will take this information from the specifications of the probe. These
    % are in um here, but the absolute scaling doesn't really matter in the
    % algorithm.

    % xcoords = 20 * [NaN NaN  1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
    % ycoords = 20 * [NaN NaN  7 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 ...
    %     17 17 18 18 19 19 20 20 21 21 22 22 23 23 24];
    xcoords = P.mrSiteXY(:,1)';
    ycoords = P.mrSiteXY(:,2)';

    % Often, multi-shank probes or tetrodes will be organized into groups of
    % channels that cannot possibly share spikes with the rest of the probe. This helps
    % the algorithm discard noisy templates shared across groups. In
    % this case, we set kcoords to indicate which group the channel belongs to.
    % In our case all channels are on the same shank in a single group so we
    % assign them all to group 1.

    % kcoords = [NaN NaN 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    if isfield(P, 'shank')
        kcoords = P.shank;
    else
        kcoords = ones(size(chanMap));
    end
    kcoords(P.viSiteZero) = nan;

    % at this point in Kilosort we do data = data(connected, :), ycoords =
    % ycoords(connected), xcoords = xcoords(connected) and kcoords =
    % kcoords(connected) and no more channel map information is needed (in particular
    % no "adjacency graphs" like in KlustaKwik).
    % Now we can save our channel map for the eMouse.

    % would be good to also save the sampling frequency here
    % fs = 25000;
    fs = P.sRateHz;

    S_chanMap = struct('chanMap', chanMap, 'connected', connected, 'xcoords', xcoords, 'ycoords', ycoords, 'kcoords', kcoords, 'fs', fs);

    struct_save_(S_chanMap, fullfile(fpath, 'chanMap.mat'));

end %func