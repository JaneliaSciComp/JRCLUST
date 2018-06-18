%--------------------------------------------------------------------------
function viTime1 = recenter_spk_(mrWav, viTime, viSite, P)
    spkLim = [-1,1] * abs(P.spkLim(1));
    viTime0 = [spkLim(1):spkLim(end)]'; %column
    miTime = bsxfun(@plus, int32(viTime0), int32(viTime(:)'));
    miTime = min(max(miTime, 1), size(mrWav, 1));
    miSite = repmat(viSite(:)', numel(viTime0), 1);
    mrWav_spk = mrWav(sub2ind(size(mrWav), miTime, miSite));
    [~, viMin_spk] = min(mrWav_spk,[],1);
    viTime_off = int32(gather_(viMin_spk') + spkLim(1) - 1);
    viTime1 = viTime + viTime_off;
    % disp(mean(viTime_off~=0))
end %func
