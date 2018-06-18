%--------------------------------------------------------------------------
function [viTime_spk, vrAmp_spk, viSite_spk] = spike_refrac_(viTime_spk, vrAmp_spk, viSite_spk, nRefrac)
    % Remove smaller spikes if a bigger one detected within nRefrac
    % spike_refrac_(viSpk, vrSpk, [], nRefrac)
    % spike_refrac_(viSpk, vrSpk, viSite, nRefrac)
    nRepeat_max = 10;
    for iRepeat = 1:nRepeat_max
        if isempty(viTime_spk), return ;end
        vnDiff_spk = diff(diff(viTime_spk) <= nRefrac);
        viStart = find(vnDiff_spk > 0) + 1;
        viEnd = find(vnDiff_spk < 0) + 1;
        if isempty(viStart) || isempty(viEnd), return; end
        viEnd(viEnd < viStart(1)) = [];
        if isempty(viEnd), return; end
        viStart(viStart > viEnd(end)) = [];
        nGroup = numel(viStart);
        if nGroup==0, return; end
        assert_(nGroup == numel(viEnd), 'spike_refrac_:nStart==nEnd');
        vlRemove = false(size(viTime_spk));

        for iGroup = 1:nGroup
            vi_ = viStart(iGroup):viEnd(iGroup);
            [~, imax_] = max(vrAmp_spk(vi_));
            vlRemove(vi_(imax_)) = true;
        end
        viTime_spk(vlRemove) = [];
        if ~isempty(vrAmp_spk), vrAmp_spk(vlRemove) = []; end
        if ~isempty(viSite_spk), viSite_spk(vlRemove) = []; end
    end
end %func
