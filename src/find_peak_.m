%--------------------------------------------------------------------------
% 9/13/17 JJJ: Edge case error fixed (vi2 indexing) for small block many loads
% 8/17/17 JJJ: fix for missing spikes
function viSpk1 = find_peak_(vrWav1, thresh1, nneigh_min)
    % nneigh_min: number of neighbors around the spike below the threshold
    %  0,1,2. # neighbors of minimum point below negative threshold
    % thresh1: absolute value. searching for negative peaks only

    if nargin<3, nneigh_min = []; end
    if isempty(nneigh_min), nneigh_min = 1; end

    viSpk1 = [];
    if isempty(vrWav1), return; end
    vl1 = vrWav1 < -abs(thresh1);
    vi2 = find(vl1);
    %vi2 = find(vrWav1 < -thresh1);
    if isempty(vi2), return; end

    if vi2(1)<=1
        if numel(vi2) == 1, return; end
        vi2(1) = [];
    end
    if vi2(end)>=numel(vrWav1)
        if numel(vi2) == 1, return; end
        vi2(end) = [];
    end
    vrWav12 = vrWav1(vi2);
    viSpk1 = vi2(vrWav12 <= vrWav1(vi2+1) & vrWav12 <= vrWav1(vi2-1));
    if isempty(viSpk1), return; end
    % viSpk1 = vi2(find(diff(diff(vrWav1(vi2))>0)>0) + 1); % only negative peak

    % if viSpk1(1) <= 1, viSpk1(1) = 2; end
    % if viSpk1(end) >= numel(vrWav1), viSpk1(end) = numel(vrWav1)-1; end
    switch nneigh_min
        case 1
        viSpk1 = viSpk1(vl1(viSpk1-1) | vl1(viSpk1+1));
        %         viSpk1 = viSpk1(vrWav1(viSpk1-1) < -thresh1 | vrWav1(viSpk1+1) < -thresh1);
        case 2
        viSpk1 = viSpk1(vl1(viSpk1-1) & vl1(viSpk1+1));
        %         viSpk1 = viSpk1(vrWav1(viSpk1-1) < -thresh1 & vrWav1(viSpk1+1) < -thresh1);
    end
end %func
