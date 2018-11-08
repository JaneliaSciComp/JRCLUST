%function [data, chanMeans, dimm_wav] = load_file_(filename, nSamples_load1, P)
function [data, chanMeans] = readRecording(filename, dtype, dshape, offset, hCfg)
    %LOADRECORDING Read a raw recording into memory as int16

    % load file to memory. output int16
    % assume that the file is chan x time
    % returns chans x time for efficient processing
    % vrWav_mean1: average across chan output

%     if ischar(filename)
%         [filename, nBytes_file1] = fopen_(vcFile, 'r');
%         if ~isempty(get_(P, 'header_offset'))
%             nBytes_file1 = nBytes_file1 - P.header_offset;
%             fseek(filename, P.header_offset, 'bof');
%         end
%         nSamples_load1 = floor(nBytes_file1 / P.nChans / bytesPerSample_(P.vcDataType));
%     else
%         vcFile = [];
%     end

    data = jrclust.utils.readBin(filename, dtype, dshape, offset);

    if strcmp(dtype, 'uint16')
        data = int16(single(data) - 2^15);
    elseif jrclust.utils.typeBytes(dtype) >= 4 % single or double
        data = int16(data/uV_per_bit);
    end

    % flip the polarity
    if hCfg.fInverse_file
        data = -data;
    end

    % extract channels
    if hCfg.fTranspose_bin % channels x samples
        if ~isempty(hCfg.chanMap)
            data = data(hCfg.chanMap, :);
        end

        chanMeans = single(mean(data, 1)); % performance improvement? (TODO: test)
        data = data';
    else % samples x channels
        if ~isempty(hCfg.chanMap)
            data = data(:, hCfg.chanMap);
        end

        % TODO: this may need pulling out of the conditional
        if ~isempty(hCfg.tlim_load) % load data only within time bounds
            nSamples = size(data, 1);

            bounds = min(max(round(hCfg.tlim_load * hCfg.sRateHz), 1), nSamples);
            data = data(bounds(1):bounds(end), :);
        end
        chanMeans = single(mean(data, 2)');
    end
end