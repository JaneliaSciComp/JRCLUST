%--------------------------------------------------------------------------
% 8/1/17 JJJ: Bugfix: If S_clu is empty, reload _jrc.mat file
function [S0, P] = load_cached_(P, fLoadWav)
    % Load cached data either from RAM or disk
    % Usage
    % S0 = load_cached_(P)
    % S0 = load_cached_(vcFile)
    % P0: Previously stored P in S0, S0.P = P overwritten
    if nargin<2, fLoadWav=1; end

    global tnWav_spk tnWav_raw trFet_spk %spike waveform (filtered)
    if ischar(P), P = loadParam_(P); end
    S0 = get(0, 'UserData');
    fClear_cache = 1;
    if ~isempty(S0)
        if isfield(S0, 'P')
            if strcmpi(P.vcFile_prm, S0.P.vcFile_prm)
                fClear_cache = 0;
            end
        end
    end
    if fClear_cache
        S0 = []; tnWav_spk = []; tnWav_raw = []; trFet_spk = []; % clear all
    end

    % Load from disk
    try
        fLoad0 = 0;
        if isempty(S0)
            fLoad0 = 1;
        elseif isempty(get_(S0, 'S_clu'))
            fLoad0 = 1;
        end

        vcFile_jrc = jrclust.utils.subsExt(P.vcFile_prm, '_jrc.mat');
        if ~exist_file_(vcFile_jrc), S0.P=P; return; end
        if fLoad0, S0 = load0_(vcFile_jrc); end
        if isempty(S0), S0.P = []; end
        [P0, S0.P] = deal(S0.P, P); %swap
        if isempty(tnWav_spk) || isempty(tnWav_raw) || isempty(trFet_spk)
            if ~fLoadWav, return; end
            if isempty(S0), return; end %no info
            try
                if get_set_(P, 'fRamCache', 1)
                    trFet_spk = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkfet.jrc'), 'single', S0.dimm_fet);
                    tnWav_spk = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkwav.jrc'), 'int16', S0.dimm_spk);
                    tnWav_raw = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkraw.jrc'), 'int16', S0.dimm_raw);
                else
                    trFet_spk = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkfet.jrc'), 'single', S0.dimm_fet);
                end
            catch
                disperr_();
            end
        end
        P = upgrade_param_(S0, P0);
        S0.P = P;
    catch hErr
        disperr_('load_cached_ error', hErr);
    end
    % S0.P = P; % force S0 to take the latest values from the prm file
    % set(0, 'UserData', S0);
end
