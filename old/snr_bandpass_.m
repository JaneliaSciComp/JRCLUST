%--------------------------------------------------------------------------
% Strip comments from cell string
% 7/24/17 JJJ: Code cleanup
function S_bp = snr_bandpass_(P)
    % compute SNR after bandpass
    % read IMEC meta file and calibrate the file. Load the first bit of the file
    % and bnad pass and estimate the SNR.
    % vrRmsQ_bp_site
    % vrRmsQ_bp_site
    % Compute SNR and apply gain correction
    % show site
    % vrRmsQ_site
    % vrSnr_clu
    % apply a correction factor

    vcFile_gain = get_set_(P, 'vcFile_gain');
    nSites = numel(P.viSite2Chan);
    if isempty(vcFile_gain)
        vr_Scale_uv_site = repmat(P.uV_per_bit, [nSites, 1]);
    else
        try
            vrGainCorr = csvread(vcFile_gain, 1, 1);
        catch
            vrGainCorr = csvread(vcFile_gain, 0, 0);
        end
    end

end %func
