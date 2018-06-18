%--------------------------------------------------------------------------
function [viSpk, vrSpk, viSite] = spike_refrac__(viSpk, vrSpk, viSite, nRefrac)
    % Remove smaller spikes if a bigger one detected within nRefrac
    % spike_refrac_(viSpk, vrSpk, [], nRefrac)
    % spike_refrac_(viSpk, vrSpk, viSite, nRefrac)

    nSkip_refrac = 8;
    % remove refractory period
    vlKeep = true(size(viSpk));
    % if isGpu_(viSpk), vlKeep = gpuArray(vlKeep); end
    while (1)
        viKeep1 = find(vlKeep);
        viRefrac1 = find(diff(viSpk(viKeep1)) <= nRefrac);
        if isempty(viRefrac1), break; end

        vi1 = viRefrac1(1:nSkip_refrac:end);
        viRemoveA = viKeep1(vi1);
        viRemoveB = viKeep1(vi1+1);
        if ~isempty(vrSpk)
            vl1 = abs(vrSpk(viRemoveA)) < abs(vrSpk(viRemoveB));
            vlKeep(viRemoveA(vl1)) = 0;
            vlKeep(viRemoveB(~vl1)) = 0;
        else
            vlKeep(viRemoveB) = 0;
        end
    end

    viSpk(~vlKeep) = [];
    if ~isempty(vrSpk), vrSpk(~vlKeep) = []; end
    if ~isempty(viSite), viSite(~vlKeep) = []; end
end %func
