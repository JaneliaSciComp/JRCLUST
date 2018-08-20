%--------------------------------------------------------------------------
function [rawTraces, channelMeans, sampleDims] = load_file_(fidBinary, nSamples, P)
    % load file to memory. output int16
    % assume that the file is chan x time
    % returns chans x time for efficient processing
    % channelMeans: average across chan output
    if ischar(fidBinary)
        filename = fidBinary;
        [fidBinary, nBytes] = fopenInfo(filename, 'r');

        if ~isempty(get_(P, 'headerOffset'))
            nBytes = nBytes - P.headerOffset;
            fseek(fidBinary, P.headerOffset, 'bof');
        end

        nSamples = floor(nBytes / (P.nChans*bytesPerSample_(P.dataType)));
    else
        filename = [];
    end

    fSingle = 0; %output single

    if P.fTranspose_bin
        sampleDims = [P.nChans, nSamples];
    else %Catalin's format
        sampleDims = [nSamples, P.nChans];
    end

    rawTraces = fread_(fidBinary, sampleDims, P.dataType);

    switch(P.dataType)
        case 'uint16'
            rawTraces = int16(single(rawTraces) - 2^15);
        case {'single', 'double'}
            rawTraces = int16(rawTraces / P.uV_per_bit);
    end

    if get_(P, 'fInverse_file')
        rawTraces = -rawTraces;
    end %flip the polarity

    % extract channels
    if P.fTranspose_bin
        if ~isempty(P.chanMap)
            rawTraces = rawTraces(P.chanMap, :);
        end

        channelMeans = single(mean(rawTraces, 1)); % 6x faster to transpose in dimm1

        if fSingle
            rawTraces = single(rawTraces') * P.uV_per_bit;
        else
            rawTraces = rawTraces';
        end
    else % Catalin's format. time x nChans
        if ~isempty(P.chanMap)
            rawTraces = rawTraces(:, P.chanMap);
        end

        if ~isempty(P.tlim_load)
            nSamples = size(rawTraces,1);
            nlim_load = min(max(round(P.tlim_load * P.sampleRateHz), 1), nSamples);
            rawTraces = rawTraces(nlim_load(1):nlim_load(end), :);
        end

        channelMeans = single(mean(rawTraces, 2)');

        if fSingle
            rawTraces = single(rawTraces) * P.uV_per_bit;
        end
    end

    if ~isempty(filename)
        fclose(fidBinary);
    end
end %func
