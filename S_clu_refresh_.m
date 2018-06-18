%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_refresh_(S_clu, fRemoveEmpty)

    if nargin<2, fRemoveEmpty=1; end
    nClu = double(max(S_clu.viClu));
    S_clu.nClu = nClu;
    viSite_spk = get0_('viSite_spk');
    % if isfield(S_clu, 'viSpk_shank'), viSite_spk = viSite_spk(S_clu.viSpk_shank); end
    % gviClu = gpuArray_(S_clu.viClu);
    % S_clu.cviSpk_clu = arrayfun(@(iClu)gather(find(gviClu==iClu)), 1:nClu, 'UniformOutput', 0);
    S_clu.cviSpk_clu = arrayfun(@(iClu)find(S_clu.viClu==iClu), 1:nClu, 'UniformOutput', 0);
    S_clu.vnSpk_clu = cellfun(@numel, S_clu.cviSpk_clu);
    S_clu.viSite_clu = double(arrayfun(@(iClu)mode(viSite_spk(S_clu.cviSpk_clu{iClu})), 1:nClu));
    if fRemoveEmpty, [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu); end
end %func
