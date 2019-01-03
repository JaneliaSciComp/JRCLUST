%--------------------------------------------------------------------------
function S0 = keyPressFcn_FigWav_(hFigWav, hEvent, S0) %amp dist
%     if nargin<3, S0 = get(0, 'UserData'); end
%     hCfg = S0.P; hClust = S0.S_clu;
%     hCfg.LineStyle=[];
    nSites = numel(hCfg.viSite2Chan);
    figData = get(hFigWav, 'UserData');

    switch lower(hEvent.Key)
        case {'uparrow', 'downarrow'}
        rescale_FigWav_(hEvent, S0, hCfg);
        clu_info_(S0); %update figpos

        case {'leftarrow', 'rightarrow', 'home', 'end'}
        % switch the current clu
        if strcmpi(hEvent.Key, 'home')
            S0.iCluCopy = 1;
        elseif strcmpi(hEvent.Key, 'end')
            S0.iCluCopy = hClust.nClu;
        elseif ~keyMod(hEvent, 'shift');
            if strcmpi(hEvent.Key, 'leftarrow')
                if S0.iCluCopy == 1, return; end
                S0.iCluCopy = S0.iCluCopy - 1;
            else
                if S0.iCluCopy == hClust.nClu, return; end
                S0.iCluCopy = S0.iCluCopy + 1;
            end
        else
            if isempty(S0.iCluPaste)
                S0.iCluPaste = S0.iCluCopy;
            end
            if strcmpi(hEvent.Key, 'leftarrow')
                if S0.iCluPaste == 1, return; end
                S0.iCluPaste = S0.iCluPaste - 1;
            else
                if S0.iCluPaste == hClust.nClu, return; end
                S0.iCluPaste = S0.iCluPaste + 1;
            end
        end
        S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0); %select first clu
        if strcmpi(hEvent.Key, 'home') || strcmpi(hEvent.Key, 'end') %'z' to recenter
            S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0);
        end
        case 'm', S0 = ui_merge_(S0); % merge clusters
        case 'space'
        % auto-select nearest cluster for black
        mrWavCor = hClust.mrWavCor;
        mrWavCor(S0.iCluCopy,S0.iCluCopy) = -inf;
        [~,S0.iCluPaste] = max(mrWavCor(:,S0.iCluCopy));
        set(0, 'UserData', S0);
        button_CluWav_simulate_([], S0.iCluPaste);
        case 's', auto_split_(1, S0);
        case 'r' %reset view
        figure_wait_(1);
        axis_([0, S0.S_clu.nClu + 1, 0, numel(hCfg.viSite2Chan) + 1]);
        figure_wait_(0);
        case {'d', 'backspace', 'delete'}, S0 = ui_delete_(S0);
        case 'z' %zoom
        iClu = S0.iCluCopy;
        iSiteClu = hClust.viSite_clu(S0.iCluCopy);
        set_axis_(hFigWav, iClu+[-1,1]*6, iSiteClu+[-1,1]*(hCfg.maxSite*2+1), [0 hClust.nClu+1], [0 nSites+1]);
        case 'c', plot_FigCorr_(S0);
        case 'v', plot_FigIsi_(S0);
        case 'a', update_spikes_(S0); clu_info_(S0);
        case 'f', clu_info_(S0);
        case 'h', msgbox_(figData.helpText, 1);
        case '0', unit_annotate_([],[], 'to_delete'); % TW
        case '1', unit_annotate_([],[], 'single'); % TW
        case '2', unit_annotate_([],[], 'multi'); % TW
        case 'numpad0', unit_annotate_([],[], 'to_delete'); % TW
        case 'numpad1', unit_annotate_([],[], 'single'); % TW
        case 'numpad2', unit_annotate_([],[], 'multi'); % TW
        case 'w', toggleVisible_(figData.hSpkAll); %toggle spike waveforms
        case 't', plot_FigTime_(S0); % time view
        case 'j', plot_FigProj_(S0); %projection view
        case 'n'
        fText = get_set_(figData, 'fText', get_set_(hCfg, 'fText', 1));
        setFigWavXTicks(figData, hClust, ~fText);
        case 'i', plot_FigHist_(S0); %ISI histogram
        case 'e', plot_FigMap_(S0);
        case 'u', update_FigCor_(S0);
        case 'p' %PSTH plot
        if isempty(hCfg.vcFile_trial), msgbox_('''vcFile_trial'' not set. Reload .prm file after setting (under "File menu")'); return; end
        plot_raster_(S0, 1);
        otherwise, figure_wait_(0); %stop waiting
    end
    figure_(hFigWav); %change the focus back to the current object
end %func
