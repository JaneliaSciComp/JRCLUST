%--------------------------------------------------------------------------
% 122917 JJJ: Got rid of Tab which is slow
function plot_raster_(S0, fNewFig)
    %plot_raster_()
    %   plot if window open using curretnly selected clusters
    %plot_raster_(P, iClu, S_clu)
    %   Open window and plot specific clusters and S_clu

    persistent hFig hFig_b
    if nargin<2, fNewFig = 0; end
    % import  trial time
    % P = loadParam(vcFile_prm);
    if ~isvalid_(hFig) && ~fNewFig, return; end
    if nargin<1, S0 = get(0, 'UserData'); end
    [P, S_clu, iCluCopy, iCluPaste] = deal(S0.P, S0.S_clu, S0.iCluCopy, S0.iCluPaste);
    if isfield(P, 'vcFile_psth'), P.vcFile_trial = P.vcFile_psth; end % old field name

    try
        if ~exist_file_(P.vcFile_trial), P.vcFile_trial = subsDir_(P.vcFile_trial, P.vcFile_prm); end
        if ~exist_file_(P.vcFile_trial)
            msgbox_(sprintf('File does not exist: vcFile_trial=%s', P.vcFile_trial), 1);
            return;
        end
        crTime_trial = loadTrial_(P.vcFile_trial);
    catch
        return;
    end
    if ~iscell(crTime_trial), crTime_trial = {crTime_trial}; end
    nstims = numel(crTime_trial);
    if isempty(crTime_trial), msgbox('Trial file does not exist', 'modal'); return; end

    [hFig, hFig_b] = create_figure_psth_(hFig, hFig_b, P, nstims);
    plot_figure_psth_(hFig, iCluCopy, crTime_trial, S_clu, P);
    if ~isempty(iCluPaste)
        set(hFig_b, 'Visible', 'on');
        plot_figure_psth_(hFig_b, iCluPaste, crTime_trial, S_clu, P);
    else
        set(hFig_b, 'Visible', 'off');
    end
end %func
