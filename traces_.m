%--------------------------------------------------------------------------
% 6/23/JJJ

%--------------------------------------------------------------------------
function traces_(P, fDebug_ui_, vcFileId, fPlot_lfp)
    % show raw traces
    % If file format is nChans x nSamples, load subset of file (fTranspose=1)
    % If file format is nSamples x nChans, load all and save to global (fTranspose=0)
    % 2017/6/22 James Jun: Added multiview (nTime_traces )
    global fDebug_ui mnWav mnWav1 % only use if P.fTranspose=0
    if nargin<4, fPlot_lfp = 0; end
    if nargin==0
        P = get0_('P');
    else
        set0_(P);
    end
    if nargin<2, fDebug_ui_=0; end
    if nargin<3, vcFileId=''; end
    if isempty(P), disperr_('traces_: P is empty'); return; end
    % S0 = load0_(P);
    fDebug_ui = fDebug_ui_;
    if ~fPlot_lfp
        S0 = load_cached_(P, 0);
        set(0, 'UserData', S0);
        set0_(fDebug_ui);
    else
        S0 = struct('P', P);
        set(0, 'UserData', S0);
    end
    % S0 = load_cached_(P, 0); %don't load raw waveform

    % get file to show
    iFile_show = 1; %files to display for clustered together
    if ~isempty(P.csFile_merge)
        csFiles_bin = filter_files_(P.csFile_merge);
        if numel(csFiles_bin)==1
            vcFile_bin = csFiles_bin{1};
        else %show multiple files
            if isempty(vcFileId)
                arrayfun(@(i)fprintf('%d: %s\n', i, csFiles_bin{i}), 1:numel(csFiles_bin), 'UniformOutput', 0);
                fprintf('---------------------------------------------\n');
                vcFileId = input('Please specify Fild ID from the list above:', 's');
            end
            if isempty(vcFileId), return; end
            iFile_show = str2num(vcFileId);
            try
                vcFile_bin = csFiles_bin{iFile_show};
            catch
                return;
            end
        end
    else
        vcFile_bin = P.vcFile; % if multiple files exist, load first
    end
    set0_(iFile_show);
    tlim_bin = P.tlim;

    % Open file
    fprintf('Opening %s\n', vcFile_bin);
    [fid_bin, nBytes_bin] = fopen_(vcFile_bin, 'r');
    if isempty(fid_bin), fprintf(2, '.bin file does not exist: %s\n', vcFile_bin); return; end
    nSamples_bin = floor(nBytes_bin / bytesPerSample_(P.vcDataType) / P.nChans);
    nLoad_bin = min(round(diff(tlim_bin) * P.sRateHz), nSamples_bin);
    if tlim_bin(1)>0
        iSample_bin = ceil(tlim_bin(1) * P.sRateHz) + 1; %offset sample number
    else
        iSample_bin = 1; %sample start location
    end
    nlim_bin = [0,nLoad_bin-1] + iSample_bin;
    if nlim_bin(1) < 1, nlim_bin = [1, nLoad_bin]; end
    if nlim_bin(2) > nSamples_bin, nlim_bin = [-nLoad_bin+1, 0] + nSamples_bin; end

    nTime_traces = get_(P, 'nTime_traces');
    [cvn_lim_bin, viRange_bin] = sample_skip_(nlim_bin, nSamples_bin, nTime_traces);
    if P.fTranspose_bin
        mnWav = [];
        fseek_(fid_bin, iSample_bin, P);
        if nTime_traces > 1
            mnWav1 = load_bin_multi_(fid_bin, cvn_lim_bin, P)';
        else
            mnWav1 = load_bin_(fid_bin, P.vcDataType, [P.nChans, nLoad_bin])'; %next keypress: update tlim_show
        end
        %     @TODO: load from cvn_lim_bin specifiers. check for end or beginning when keyboard command
    else %load whole thing
        mnWav = load_bin_(fid_bin, P.vcDataType, [nSamples_bin, P.nChans]); %next keypress: update tlim_show
        fclose(fid_bin);
        fid_bin = [];
        %mnWav1 = mnWav((nlim_bin(1):nlim_bin(2)), :);
        mnWav1 = mnWav(viRange_bin, :);
        disp('Entire raw traces are cached to RAM since fTranspose=0.');
    end %if
    mnWav1 = uint2int_(mnWav1);

    % full screen width
    hFig_traces = create_figure_('Fig_traces', [0 0 .5 1], vcFile_bin, 0, 1); %remove all other figure traces
    hAx = axes_new_(hFig_traces); % create axis
    hPlot = line(hAx, nan, nan, 'Color', [1 1 1]*.5, 'Parent', hAx, 'LineWidth', .5);
    hPlot_edges = plot(nan, nan, 'Color', [1 0 0]*.5, 'Parent', hAx, 'LineWidth', 1);
    set(hAx, 'Position',[.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual');
    S_fig = makeStruct_(hAx, hPlot, nlim_bin, fid_bin, nSamples_bin, nLoad_bin, hPlot_edges);
    S_fig.maxAmp = P.maxAmp;
    S_fig.vcTitle = '[H]elp; (Sft)[Up/Down]:Scale(%0.1f uV); (Sft)[Left/Right]:Time; [F]ilter; [J]ump T; [C]han. query; [R]eset view; [P]SD; [S]pike; [A]ux chan; [E]xport; [T]race; [G]rid';
    S_fig.csHelp = { ...
    'Left/Right: change time (Shift: x4)', ...
    '[J]ump T', ...
    '[Home/End]: go to beginning/end of file', ...
    '---------', ...
    'Up/Down: change scale (Shift: x4)', ...
    'Zoom: Mouse wheel', ...
    '[x/y/ESC]: zoom direction', ...
    'Pan: hold down the wheel and drag', ...
    '[R]eset view', ...
    '---------', ...
    '[F]ilter toggle', ...
    '[S]pike toggle', ...
    'Gri[D] toggle', ...
    '[T]races toggle', ...
    '---------', ...
    '[C]hannel query', ...
    '[A]ux channel display', ...
    '[P]ower spectrum', ...
    '[E]xport to workspace', ...
    };
    S_fig = struct_append_(S_fig, ...
    struct('vcGrid', 'on', 'vcFilter', 'off', 'vcSpikes', 'on', 'vcTraces', 'on'));
    set(hFig_traces, 'UserData', S_fig);
    set(hFig_traces, 'color', 'w', 'KeyPressFcn', @keyPressFcn_Fig_traces_, 'BusyAction', 'cancel', 'CloseRequestFcn', @close_hFig_traces_);
    mouse_figure(hFig_traces);

    Fig_traces_plot_(1); % Plot spikes and color clusters
end %func
