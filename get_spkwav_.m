%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function tnWav_ = get_spkwav_(P, fRaw)
    % if ~fRamCache, only keep one of the tnWav_raw or tnWav_spk in memory
    global tnWav_spk tnWav_raw
    if nargin<1, P = []; end
    if isempty(P), P = get0_('P'); end
    if nargin<2, fRaw = P.fWav_raw_show; end

    fRamCache = get_set_(P, 'fRamCache', 1);
    if fRaw
        if ~fRamCache, tnWav_spk = []; end % clear spk
        if isempty(tnWav_raw), tnWav_raw = load_spkraw_(); end
        tnWav_ = tnWav_raw;
    else
        if ~fRamCache, tnWav_raw = []; end % clear raw
        if isempty(tnWav_spk), tnWav_spk = load_spkwav_(); end
        tnWav_ = tnWav_spk;
    end
end %func
