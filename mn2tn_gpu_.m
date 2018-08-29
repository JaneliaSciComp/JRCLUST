%--------------------------------------------------------------------------
function [tr, miRange] = mn2tn_gpu_(mr, spkLim, viTime, viSite)
    % what to do if viTime goves out of the range?
    % gpu memory efficient implementation
    % it uses GPU if mr is in GPU

    if nargin<4, viSite=[]; end %faster indexing
    % if nargin<5, fMeanSubt=0; end

    % JJJ 2015 Dec 24
    % vr2mr2: quick version and doesn't kill index out of range
    % assumes vi is within range and tolerates spkLim part of being outside
    % works for any datatype
    if isempty(viTime), tr=[]; return; end
    [nT, nSites] = size(mr);
    if ~isempty(viSite)
        mr = mr(:, viSite);
        nSites = numel(viSite);
    end
    if iscolumn(viTime), viTime = viTime'; end

    useGPU = isGpu_(mr);
    viTime = gpuArray_(viTime, useGPU);
    spkLim = gpuArray_(spkLim, useGPU);

    viTime0 = [spkLim(1):spkLim(end)]'; %column
    miRange = bsxfun(@plus, int32(viTime0), int32(viTime));
    miRange = min(max(miRange, 1), nT);
    % miRange = miRange(:);
    tr = zeros([numel(viTime0), numel(viTime), nSites], 'int16');
    dimm_tr = size(tr);
    for iSite = 1:nSites
        if useGPU
            %         vr1 = gpuArray(mr(:,iSite));
            %         tr(:,:,iSite) = gather(vr1(miRange));
            tr(:,:,iSite) = gather_(reshape(mr(miRange, iSite), dimm_tr(1:2)));
        else
            tr(:,:,iSite) = reshape(mr(miRange, iSite), dimm_tr(1:2));
        end
    end
end % function
