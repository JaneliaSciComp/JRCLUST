%--------------------------------------------------------------------------
% 6/23/JJJ

%--------------------------------------------------------------------------
function traces_(hCfg, recID, showLFP)
    % show raw traces
    % If file format is nChans x nSamples, load subset of file (fTranspose=1)
    % If file format is nSamples x nChans, load all and save to global (fTranspose=0)
    % 2017/6/22 James Jun: Added multiview (nTime_traces )
    global mnWav mnWav1 % only use if P.fTranspose=0
    if nargin < 4
        showLFP = 0;
    end
    if nargin < 3
        recID = '';
    end

    if ~showLFP
        S0 = load_cached_(hCfg, 0);
        set(0, 'UserData', S0);
    else
        S0 = struct('P', hCfg);
        set(0, 'UserData', S0);
    end

    % get file to show
    iFile_show = 1; %files to display for clustered together
    if numel(hCfg.rawRecordings) > 1
        if isempty(recID)
            arrayfun(@(i) fprintf('%d: %s\n', i, hCfg.rawRecordings{i}), 1:numel(hCfg.rawRecordings), 'UniformOutput', 0);
            fprintf('---------------------------------------------\n');
            recID = input('Please specify File ID from the list above:', 's');
        end

        if isempty(recID)
            return;
        end

        iFile_show = str2double(recID);

        try
            recFilename = hCfg.rawRecordings{iFile_show};
        catch
            return;
        end
    else
        recFilename = hCfg.rawRecordings{1}; % if multiple files exist, load first
    end

    set0_(iFile_show);

    % Open file
    fprintf('Opening %s\n', recFilename);
    hRec = jrclust.models.recording.Recording(recFilename, hCfg.dtype, hCfg.nChans, hCfg.headerOffset, hCfg);

%     [fid_bin, nBytes_bin] = fopen_(vcFile_bin, 'r');
%     if isempty(fid_bin)
%         fprintf(2, '.bin file does not exist: %s\n', vcFile_bin);
%         return;
%     end
%     nSamples_bin = floor(nBytes_bin / bytesPerSample_(hCfg.vcDataType) / hCfg.nChans);
    nSamplesTotal = hRec.nSamples;
    windowWidth = min(floor(diff(hCfg.dispTimeLimits) * hCfg.sampleRate), nSamplesTotal);

    if hCfg.dispTimeLimits(1) > 0
        iSampleOffset = ceil(hCfg.dispTimeLimits(1) * hCfg.sampleRate) + 1; % offset sample number
    else
        iSampleOffset = 1; %sample start location
    end

    windowBounds = [0, windowWidth - 1] + iSampleOffset;
    if windowBounds(1) < 1
        windowBounds = [1, windowWidth];
    end
    if windowBounds(2) > nSamplesTotal
        windowBounds = [-windowWidth + 1, 0] + nSamplesTotal;
    end

    [cvn_lim_bin, viRange_bin] = sample_skip_(windowBounds, nSamplesTotal, hCfg.nTime_traces);

    mnWav = [];
    mnWav1 = cellfun(@(lims) hRec.readROI(hCfg.siteMap, lims(1):lims(2)), cvn_lim_bin, 'UniformOutput', false);
    mnWav1 = jrclust.utils.neCell2mat(mnWav1);

%     if hCfg.fTranspose_bin
%         mnWav = [];
%         fseek_(fid_bin, iSampleOffset, hCfg);
%         if hCfg.nTime_traces > 1
%             mnWav1 = load_bin_multi_(fid_bin, cvn_lim_bin, hCfg)';
%         else
%             mnWav1 = jrclust.utils.readBin(fid_bin, hCfg.vcDataType, [hCfg.nChans, windowWidth])'; %next keypress: update tlim_show
%         end
%         %     @TODO: load from cvn_lim_bin specifiers. check for end or beginning when keyboard command
%     else %load whole thing
%         mnWav = jrclust.utils.readBin(fid_bin, hCfg.vcDataType, [nSamplesTotal, hCfg.nChans]); %next keypress: update tlim_show
%         fclose(fid_bin);
%         fid_bin = [];
%         %mnWav1 = mnWav((nlim_bin(1):nlim_bin(2)), :);
%         mnWav1 = mnWav(viRange_bin, :);
%         disp('Entire raw traces are cached to RAM since fTranspose=0.');
%     end %if
    mnWav1 = uint2int_(mnWav1);

    % full screen width
    hFigTraces = jrclust.views.Figure('FigTraces', [0 0 .5 1], recFilename, false, true);
    %hFigTraces = create_figure_('Fig_traces', [0 0 .5 1], recFilename, 0, 1); %remove all other figure traces
    hFigTraces.axes(); % create axis
    %hPlot = line(hAx, nan, nan, 'Color', [1 1 1]*.5, 'Parent', hAx, 'LineWidth', .5);
    hFigTraces.addPlot('hLine', @line, nan, nan, 'Color', [1 1 1]*.5, 'LineWidth', .5);
    %hPlot_edges = plot(nan, nan, 'Color', [1 0 0]*.5, 'Parent', hAx, 'LineWidth', 1);
    hFigTraces.addPlot('hEdges', nan, nan, 'Color', [1 0 0]*.5, 'LineWidth', 1);
    %set(hAx, 'Position',[.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual');
    hFigTraces.axApply(@set, 'Position', [.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual');
    %figData = makeStruct_(hAx, hPlot, windowBounds, fid_bin, nSamplesTotal, windowWidth, hPlot_edges);
    figData.maxAmp = hCfg.maxAmp;
    figData.vcTitle = '[H]elp; (Sft)[Up/Down]:Scale(%0.1f uV); (Sft)[Left/Right]:Time; [F]ilter; [J]ump T; [C]han. query; [R]eset view; [P]SD; [S]pike; [A]ux chan; [E]xport; [T]race; [G]rid';
    figData.helpText = { ...
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
    figData = struct_append_(figData, ...
        struct('vcGrid', 'on', 'vcFilter', 'off', 'vcSpikes', 'on', 'vcTraces', 'on'));
    hFigTraces.figData = figData;
    hFigTraces.figSet('color', 'w', 'BusyAction', 'cancel', 'CloseRequestFcn', @close_hFig_traces_);
    hFigTraces.keyPressFcn = @keyPressFcn_Fig_traces_;
    hFigTraces.setMouseable();

    doPlotFigTraces(1); % Plot spikes and color clusters
end
