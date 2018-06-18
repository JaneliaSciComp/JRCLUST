%--------------------------------------------------------------------------
% 10/10/17 JJJ: created and tested.
function tnWav_spk1 = fread_spkwav_(viSpk1, fRamCache)
    % Read spikes from file
    persistent fid_ dimm_
    if nargin<2, P = get0_('P'); fRamCache = get_set_(P, 'fRamCache', 1); end

    % reset file
    if nargin==0
        if ~isempty(fid_), fclose_(fid_, 1); end
        [fid_, tnWav_spk1] = deal([]);
        return;
    end

    if fRamCache
        global tnWav_spk
        if isempty(tnWav_spk), tnWav_spk = load_spkwav_(); end
        tnWav_spk1 = tnWav_spk(:,:,viSpk1);
        return;
    end

    % open file
    if isempty(fid_)
        [P, dimm_] = get0_('P', 'dimm_spk');
        vcFile_ = strrep(P.vcFile_prm, '.prm', '_spkwav.jrc');
        assert_(exist_file_(vcFile_), sprintf('File must exist: %s\n', vcFile_));
        fid_ = fopen(vcFile_, 'r');
    end

    tnWav_spk1 = fread_spk_(fid_, dimm_, viSpk1);
end %func
