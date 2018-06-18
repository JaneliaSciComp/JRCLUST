%--------------------------------------------------------------------------
function [vrY, vrX] = tr2plot_(trWav, iClu, viSite_show, maxAmp, P)
    if nargin<2, iClu=1; end
    iClu = double(iClu);
    if nargin<5, P = get0_('P'); end %S0 = get(0, 'UserData'); P = S0.P;
    % [~, S_fig] = get_fig_cache_('FigWav');
    % if isfield(S_fig, 'maxAmp')
    %     maxAmp = S_fig.maxAmp;
    % else
    %     maxAmp = P.maxAmp;
    % end

    if nargin<3, viSite_show = []; end
    % P = funcDefStr_(P, 'LineStyle', 'k', 'spkLim', [-10 24], 'maxAmp', 500, 'viSite_show', []);
    % P.LineStyle
    % if isempty(P.LineStyle), P.LineStyle='k'; end
    if isempty(viSite_show), viSite_show = 1:size(trWav,2); end

    [nSamples, nChans, nSpk] = size(trWav);
    nSites_show = numel(viSite_show);
    trWav = single(trWav) / maxAmp;
    trWav = trWav + repmat(single(viSite_show(:)'), [size(trWav,1),1,size(trWav,3)]);
    trWav([1,end],:,:) = nan;
    vrY = trWav(:);

    if nargout>=2
        vrX = wav_clu_x_(iClu, P);
        vrX = repmat(vrX(:), [1, nSites_show * nSpk]);
        vrX = single(vrX(:));
    end
end
