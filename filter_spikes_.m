%--------------------------------------------------------------------------
function [viTime_spk11, viSite_spk11] = filter_spikes_(viTime_spk0, viSite_spk0, tlim)
    % Filter spikes that is within tlim specified

    [viTime_spk11, viSite_spk11] = deal([]);
    if isempty(viTime_spk0), return; end
    viKeep11 = find(viTime_spk0 >= tlim(1) & viTime_spk0 <= tlim(end));
    viTime_spk11 = viTime_spk0(viKeep11)  + (1 - tlim(1)); % shift spike timing
    if ~isempty(viSite_spk0)
        viSite_spk11 = viSite_spk0(viKeep11);
    else
        viSite_spk11 = [];
    end
end %func
