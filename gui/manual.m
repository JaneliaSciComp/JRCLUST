%--------------------------------------------------------------------------
function manual(P, debugMode, alg)
    % display manual sorting interface
    global fDebug_ui trFet_spk

    if nargin < 2
        debugMode = 0;
    end
    
    if nargin < 3
        alg = 'JRCLUST';
    end

    % Load info
    if ~isSorted(P)
        error(['File must to be sorted first (run "jrc spikesort "', P.prmFile, '")']);
    end

    [S0, P] = load_cached_(P);

    if ~isfield(S0, 'mrPos_spk')
        S0.mrPos_spk = spk_pos_(S0, trFet_spk);
        set(0, 'UserData', S0);
    end

    fDebug_ui = 0;
    P.useGPU = 0; % do not use GPU for manual use
    set0_(fDebug_ui, P);

    if debugMode
        fDebug_ui = 1;
        S0 = set0_(fDebug_ui);
%         [S_clu, S0] = post_merge_(S0.S_clu, P); % redo the clustering (reset to auto)
        S0 = set0_(P);
    else
        if ~isempty(get_set_(S0, 'cS_log', {}))
            switch lower(questdlg_('Load last saved?', 'Confirmation'))
                case 'no'
                    [S_clu, S0] = post_merge_(S0.S_clu, P);
                    S0 = clear_log_(S0);
                case 'cancel'
                    return;
                case 'yes'
                    S0 = set0_(P); % update the P structure
                    S0.S_clu = S_clu_update_wav_(S0.S_clu, P);
            end
        end
    end

    % Create figures
    hMsg = msgbox_('Plotting... (this closes automatically)');
    t1 = tic;

    set(0, 'UserData', S0);
    S0 = figures_manual_(P); % create figures for manual interface

    clear mouse_figure;
    clear get_fig_cache_ get_tag_ %clear persistent figure handles

    % Set fields
    S0 = struct_merge_(S0, ...
        struct('iCluCopy', 1, 'iCluPaste', [], 'hCopy', [], 'hPaste', [], 'nSites', numel(P.viSite2Chan)));
    set(0, 'UserData', S0);

    % hFigRD
    S0.S_clu = plot_FigRD_(S0.S_clu, P); % ask user before doing so

    % Set initial amplitudes
    set(0, 'UserData', S0);
    plot_FigWavCor_(S0); % hFigWavCor
    S0 = plot_FigWav_(S0); % hFigWav %do this after for ordering

    % hFigProj, hFigHist, hFigIsi, hFigCorr, hFigPos, hFigMap, hFigTime
    tryClose(get_fig_('FigTrial')); %close previous FigTrial figure
    tryClose(get_fig_('FigTrial_b')); %close previous FigTrial figure
    S0 = button_CluWav_simulate_(1, [], S0, 1); %select first clu TW
    auto_scale_proj_time_(S0);
    S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
    %S0.cS_log = load_(strrep(P.prmFile, '.prm', '_log.mat'), 'cS_log', 0);
    S_log = load_(strrep(P.prmFile, '.prm', '_log.mat'), [], 0);
    if ~isempty(S_log), S0.cS_log = {S_log}; end
    save_log_('start', S0); %crash proof log

    % Finish up
    tryClose(hMsg);
    fprintf('UI creation took %0.1fs\n', toc(t1));
end %func
