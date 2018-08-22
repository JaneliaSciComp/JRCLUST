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
    if ~tryIsValid(hFig) && ~fNewFig, return; end
    if nargin<1, S0 = get(0, 'UserData'); end
    [P, S_clu, primarySelectedCluster, secondarySelectedCluster] = deal(S0.P, S0.S_clu, S0.primarySelectedCluster, S0.secondarySelectedCluster);
    if isfield(P, 'vcFile_psth'), P.trialFile = P.vcFile_psth; end % old field name

    try
        % begin TW block
        if ~isfield(P, 'trialFile') || isempty(P.trialFile)
            if exist(strrep(P.paramFile,".prm",".starts.mat"),'file')
                P.trialFile=char(strrep(P.paramFile,".prm",".starts.mat"));
            elseif exist(strrep(P.paramFile,".prm",".mat"),'file')
                P.trialFile=char(strrep(P.paramFile,".prm",".mat"));
            else
                msgbox_('''trialFile'' not set. Reload .prm file after setting (under "File menu")'); return;
            end
        end
        % end TW block

        if ~fileExists(P.trialFile), P.trialFile = replaceDir(P.trialFile, P.paramFile); end
        if ~fileExists(P.trialFile)
            msgbox_(sprintf('File does not exist: trialFile=%s', P.trialFile), 1);
            return;
        end
        crTime_trial = loadTrial_(P.trialFile);
    catch
        return;
    end
    if ~iscell(crTime_trial), crTime_trial = {crTime_trial}; end
    nstims = numel(crTime_trial);
    if isempty(crTime_trial), msgbox('Trial file does not exist', 'modal'); return; end

    [hFig, hFig_b] = createFigPSTH(hFig, hFig_b, P, nstims);
    plot_figure_psth_(hFig, primarySelectedCluster, crTime_trial, S_clu, P);
    if ~isempty(secondarySelectedCluster)
        set(hFig_b, 'Visible', 'on');
        plot_figure_psth_(hFig_b, secondarySelectedCluster, crTime_trial, S_clu, P);
    else
        set(hFig_b, 'Visible', 'off');
    end
end %func
