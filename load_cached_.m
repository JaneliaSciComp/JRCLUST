%--------------------------------------------------------------------------
% 8/1/17 JJJ: Bugfix: If S_clu is empty, reload _jrc.mat file
function [S0, P] = load_cached_(P, loadWaveforms)
    % Load cached data either from RAM or disk
    % Usage
    % S0 = load_cached_(P)
    % S0 = load_cached_(vcFile)
    % P0: Previously stored P in S0, S0.P = P overwritten
    if nargin < 2
        loadWaveforms = 1;
    end

    global spikeWaveforms spikeTraces spikeFeatures %spike waveform (filtered)
    if ischar(P), P = loadParams(P); end
    S0 = get(0, 'UserData');
    fClear_cache = 1;
    if ~isempty(S0)
        if isfield(S0, 'P')
            if strcmpi(P.paramFile, S0.P.paramFile)
                fClear_cache = 0;
            end
        end
    end
    if fClear_cache
        S0 = []; spikeWaveforms = []; spikeTraces = []; spikeFeatures = []; % clear all
    end

    % Load from disk
    try
        fLoad0 = 0;
        if isempty(S0)
            fLoad0 = 1;
        elseif isempty(get_(S0, 'S_clu'))
            fLoad0 = 1;
        end

        vcFile_jrc = strrep(P.paramFile, '.prm', '_jrc.mat');
        if ~fileExists(vcFile_jrc), S0.P=P; return; end
        if fLoad0, S0 = load0_(vcFile_jrc); end
        if isempty(S0), S0.P = []; end
        [P0, S0.P] = deal(S0.P, P); %swap
        if isempty(spikeWaveforms) || isempty(spikeTraces) || isempty(spikeFeatures)
            if ~loadWaveforms, return; end
            if isempty(S0), return; end %no info
            try
                if get_set_(P, 'fRamCache', 1)
                    spikeFeatures = load_bin_(strrep(P.paramFile, '.prm', '_spkfet.jrc'), 'single', S0.featureDims);
                    spikeWaveforms = load_bin_(strrep(P.paramFile, '.prm', '_spkwav.jrc'), 'int16', S0.waveformDims);
                    spikeTraces = load_bin_(strrep(P.paramFile, '.prm', '_spkraw.jrc'), 'int16', S0.traceDims);
                else
                    spikeFeatures = load_bin_(strrep(P.paramFile, '.prm', '_spkfet.jrc'), 'single', S0.featureDims);
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
end %func
