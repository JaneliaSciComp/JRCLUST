%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function plot_aux_rate_(fSelectedUnit)
    % Aux channel vs. rate
    if nargin<1, fSelectedUnit = 0; end %plot all
    [P, S_clu, iCluCopy] = get0_('P', 'S_clu', 'iCluCopy');
    P = loadParams(P.prmFile);
    [vrWav_aux, vrTime_aux] = load_aux_(P);
    if isempty(vrWav_aux), msgbox_('Aux input is not found'); return; end
    mrRate_clu = clu_rate_(S_clu, [], numel(vrWav_aux));
    vrCorr_aux_clu = arrayfun(@(i)corr(vrWav_aux, mrRate_clu(:,i), 'type', 'Pearson'), 1:size(mrRate_clu,2));
    if ~fSelectedUnit, iCluCopy = []; end
    plot_aux_corr_(mrRate_clu, vrWav_aux, vrCorr_aux_clu, vrTime_aux, iCluCopy);
    vcMsg = assignWorkspace_(mrRate_clu, vrWav_aux, vrCorr_aux_clu, vrTime_aux);
    % msgbox_(vcMsg);
end %func
