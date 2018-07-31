%--------------------------------------------------------------------------
% Make trial function migrated from jrc v1 (jrclust.m)
function make_trial_(vcFile_prm, fImec)
    % fImec: code-based event detection
    % make a _trial.mat file. stores real time
    if nargin<2, fImec = 0; end

    P = loadParams(vcFile_prm);

    % ask which channel and which output file
    csAns = inputdlg('Which channel to load', 'Channel', 1, {num2str(P.nChans)});
    iChan = str2double(csAns{1});

    % get output file
    vcFile_trial = subsFileExt(P.prmFile,  sprintf('_ch%d_trial.mat', iChan));
    try
        [FileName,PathName,FilterIndex] = uiputfile(vcFile_trial, 'Save file name');
        if ~FilterIndex, return; end %cancelled
        vcFile_trial = [PathName, FileName];
    catch
        fprintf('uiputfile error (old Matlab version). Accepting default');
    end

    vcAns = questdlg_('TTL Edge', 'Select rising or falling edge for the TTL pulses', 'Rising edge', 'Falling edge', 'Rising edge');
    if isempty(vcAns), return; end

    hMsg = msgbox_('Loading... (this closes automatically)');
    vrWav = load_bin_chan_(P, iChan);
    if isempty(vrWav), fprintf(2, 'File loading error: %s\n', P.prmFile); return; end
    % fid = memmapfile(P.vcFile, 'Offset', 0, 'Format', P.vcDataType, 'Repeat', inf);
    % vrWav = fid.Data(iChan:P.nChans:end);
    % clear fid;

    if fImec % convert to TTL signal
        dinput_imec = get_set_(P, 'dinput_imec_trial', 1);
        fprintf('Digital input %d selected (P.dinput_imec_trial: 1-16)\n', dinput_imec);
        codeval = int16(bitshift(1,dinput_imec-1));
        vrWav = single(vrWav == codeval);
    end

    [maxV, minV] = deal(max(vrWav), min(vrWav));
    if isempty(get_(P, 'thresh_trial'))
        thresh = (maxV + minV)/2;
    else
        thresh = P.thresh_trial; % load a specific threshold only
    end
    nRefrac_trial = round(P.tRefrac_trial * P.sRateHz);
    viT_rising = remove_refrac(find(vrWav(1:end-1) < thresh & vrWav(2:end) >= thresh) + 1, nRefrac_trial);
    viT_falling = remove_refrac(find(vrWav(1:end-1) > thresh & vrWav(2:end) <= thresh), nRefrac_trial);
    viT = ifeq_(strcmpi(vcAns, 'Rising edge'), viT_rising, viT_falling);
    vrT = viT / P.sRateHz;

    % save
    save(vcFile_trial, 'vrT');
    disp(['Saved to ', vcFile_trial]);
    vcTitle = sprintf('%d events detected (%s)\n', numel(vrT), lower(vcAns));

    % edit .prm file
    P_trial.vcFile_trial = vcFile_trial;
    edit_prm_file_(P_trial, vcFile_prm);
    fprintf('Parameter file updated: %s\n\tvcFile_trial = ''%s''\n', vcFile_prm, vcFile_trial);

    % plot
    hFig = createFigure('Trial timing', [0 0 .5 1], vcFile_trial, 1, 1); hold on;
    % vlOver = vr_set_(vrWav >= thresh, [viT_rising(:) - 1; viT_falling(:) + 1], 1);
    vlOver = vrWav >= thresh;
    plot(find(vlOver)/P.sRateHz, vrWav(vlOver), 'b.');
    stem(vrT, vrWav(viT), 'r'); hold on;
    plot(get(gca, 'XLim'), double(thresh) * [1 1], 'm-');
    xylabel_(gca, 'Time (s)', sprintf('Chan %d', iChan), vcTitle);
    ylim([minV, maxV]); grid on;
    tryClose(hMsg);
end
