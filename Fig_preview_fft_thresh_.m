%--------------------------------------------------------------------------
function Fig_preview_fft_thresh_(hFig)
    S_fig = get(hFig, 'UserData');
    P = get0_('P');

    % ask user
    if isempty(S_fig.fft_thresh)
        vc_ = '0';
    else
        vc_ = num2str(S_fig.fft_thresh);
    end
    vcAns = inputdlg_('Set a threshold for FFT cleanup (0 to disable)', 'fft_thresh (20 recommended)', 1, {vc_});
    if isempty(vcAns), return; end
    fft_thresh = str2double(vcAns{1});
    if isnan(fft_thresh), return; end

    S_fig.fft_thresh = fft_thresh;
    S_fig = Fig_preview_update_(hFig, S_fig, 1);
end %func
