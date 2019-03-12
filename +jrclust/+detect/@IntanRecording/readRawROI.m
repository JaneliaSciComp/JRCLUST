function roi = readRawROI(obj, rows, cols)
    %READRAWROI Get a region of interest by rows/cols from the raw file
    % How many data blocks remain in this file?
    bytesRemaining = obj.fSizeBytes - obj.headerOffset;
    dataPresent = (bytesRemaining > 0);

    if (dataPresent)
        if ~obj.rawIsOpen
            obj.openRaw();
            closeAfter = 1;
        else
            closeAfter = 0;
        end
        % Pre-allocate memory for data
        roi = zeros(numel(rows), numel(cols), 'single');
        blocks = floor((cols-1)/obj.nSamplesBlock) + 1;
        uniqueBlocks = unique(blocks); % these will be sorted
        blockOffsets = mod(cols, obj.nSamplesBlock);
        blockOffsets(blockOffsets == 0) = obj.nSamplesBlock;

        % go to just past the header
        for iiBlock = 1:numel(uniqueBlocks)
            iBlock = uniqueBlocks(iiBlock);
            blockMask = (blocks == iBlock);

            % seek to just past timestamps for this block
            fseek(obj.rawFid, obj.headerOffset + (iBlock-1)*obj.nBytesBlock + 4*obj.nSamplesBlock, 'bof');
            % read in all channels/samples for this block and subset
            blockROI = fread(obj.rawFid, [obj.nSamplesBlock, obj.nChans], 'uint16')';
            roi(:, blockMask) = 0.195 * (blockROI(rows, blockOffsets(blockMask)) - 32768);
        end

        if (obj.notchFilterFreq > 0)
            for iChan = 1:obj.nChans
                roi(iChan, :) =  notchFilter(roi(iChan, :), obj.sampleRate, obj.notchFilterFreq, 10);
            end
        end

        if closeAfter
            obj.closeRaw();
        end
    else
        roi = [];
    end
end

function out = notchFilter(in, fSample, fNotch, Bandwidth)

    % out = notch_filter(in, fSample, fNotch, Bandwidth)
    %
    % Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
    % fSample = sample rate of data (in Hz or Samples/sec)
    % fNotch = filter notch frequency (in Hz)
    % Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
    %   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
    %   poor time-domain properties with an extended ringing response to
    %   transient disturbances.
    %
    % Example:  If neural data was sampled at 30 kSamples/sec
    % and you wish to implement a 60 Hz notch filter:
    %
    % out = notch_filter(in, 30000, 60, 10);

    tstep = 1/fSample;
    Fc = fNotch*tstep;

    L = length(in);

    % Calculate IIR filter parameters
    d = exp(-2*pi*(Bandwidth/2)*tstep);
    b = (1 + d*d)*cos(2*pi*Fc);
    a0 = 1;
    a1 = -b;
    a2 = d*d;
    a = (1 + d*d)/2;
    b0 = 1;
    b1 = -2*cos(2*pi*Fc);
    b2 = 1;

    out = zeros(size(in));
    out(1) = in(1);  
    out(2) = in(2);
    % (If filtering a continuous data stream, change out(1) and out(2) to the
    %  previous final two values of out.)

    % Run filter
    for i=3:L
        out(i) = (a*b2*in(i-2) + a*b1*in(i-1) + a*b0*in(i) - a2*out(i-2) - a1*out(i-1))/a0;
    end

end