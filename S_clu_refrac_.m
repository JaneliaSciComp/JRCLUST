%--------------------------------------------------------------------------
function [S_clu, nRemoved] = S_clu_refrac_(S_clu, P, iClu1)
    % clu_refrac(Sclu, P)   %process refrac on all clusters
    % clu_refrac(Sclu, P, iClu1) %process on specific clusters
    % P.nSkip_refrac = 4;
    % P.fShow_refrac = 0;
    spikeTimes = get0_('spikeTimes');
    % remove refractory spikes
    if nargin==2
        %     P = varargin{1}; %second input
        nClu = max(S_clu.spikeClusters);
        P.fShow_refrac = 1;
        nRemoved = 0;
        for iClu=1:nClu
            [S_clu, nRemoved1] = S_clu_refrac_(S_clu, P, iClu);
            nRemoved = nRemoved + nRemoved1;
        end
        return;
    else
        %     iClu1 = varargin{1};
        %     P = varargin{2};
        nRemoved = 0;
        if ~isfield(P, 'nSkip_refrac'), P.nSkip_refrac = 4; end
        if ~isfield(P, 'fShow_refrac'), P.fShow_refrac = 1; end
        try
            viSpk1 = S_clu.cviSpk_clu{iClu1};
        catch
            viSpk1 = find(S_clu.spikeClusters == iClu1);
        end
        if isempty(viSpk1), return; end

        viTime1 = spikeTimes(viSpk1);
        nRefrac = round(P.spkRefrac_ms * P.sRateHz / 1000);

        % removal loop
        vlKeep1 = true(size(viTime1));
        while (1)
            viKeep1 = find(vlKeep1);
            viRefrac11 = find(diff(viTime1(viKeep1)) < nRefrac) + 1;
            if isempty(viRefrac11), break; end

            vlKeep1(viKeep1(viRefrac11(1:P.nSkip_refrac:end))) = 0;
        end
        nRemoved = sum(~vlKeep1);
        nTotal1 = numel(vlKeep1);
        S_clu.spikeClusters(viSpk1(~vlKeep1)) = 0;

        S_clu.cviSpk_clu{iClu1} = viSpk1(vlKeep1);
        S_clu.vnSpk_clu(iClu1) = sum(vlKeep1);
    end

    if get_(P, 'fVerbose')
        fprintf('Clu%d removed %d/%d (%0.1f%%) duplicate spikes\n', ...
        iClu1, nRemoved, nTotal1, nRemoved/nTotal1*100);
    end
end %func
