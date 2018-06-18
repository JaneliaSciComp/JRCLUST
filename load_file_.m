%--------------------------------------------------------------------------
function [mnWav1, vrWav_mean1, dimm_wav] = load_file_(fid_bin, nSamples_load1, P)
    % load file to memory. output int16
    % assume that the file is chan x time
    % returns chans x time for efficient processing
    % vrWav_mean1: average across chan output
    if ischar(fid_bin)
        vcFile = fid_bin;
        [fid_bin, nBytes_file1] = fopen_(vcFile, 'r');
        if ~isempty(get_(P, 'header_offset'))
            nBytes_file1 = nBytes_file1 - P.header_offset;
            fseek(fid_bin, P.header_offset, 'bof');
        end
        nSamples_load1 = floor(nBytes_file1 / P.nChans / bytesPerSample_(P.vcDataType));
    else
        vcFile = [];
    end
    fSingle = 0; %output single
    if P.fTranspose_bin
        dimm_wav = [P.nChans, nSamples_load1];
    else %Catalin's format
        dimm_wav = [nSamples_load1, P.nChans];
    end
    mnWav1 = fread_(fid_bin, dimm_wav, P.vcDataType);
    % [mnWav1, P.fGpu] = gpuArray_(mnWav1, P.fGpu);
    switch(P.vcDataType)
        case 'uint16', mnWav1 = int16(single(mnWav1)-2^15);
        case {'single', 'double'}, mnWav1 = int16(mnWav1 / P.uV_per_bit);
    end
    if get_(P, 'fInverse_file'), mnWav1 = -mnWav1; end %flip the polarity

    % extract channels
    if P.fTranspose_bin
        if ~isempty(P.viSite2Chan), mnWav1 = mnWav1(P.viSite2Chan,:); end
        vrWav_mean1 = single(mean(mnWav1, 1)); %6x faster to transpose in dimm1
        if fSingle
            mnWav1 = single(mnWav1') * P.uV_per_bit;
        else
            mnWav1 = mnWav1';
        end
    else %Catalin's format. time x nChans
        if ~isempty(P.viSite2Chan), mnWav1 = mnWav1(:,P.viSite2Chan); end
        if ~isempty(P.tlim_load)
            nSamples = size(mnWav1,1);
            nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
            mnWav1 = mnWav1(nlim_load(1):nlim_load(end), :);
        end
        vrWav_mean1 = single(mean(mnWav1, 2)');
        if fSingle
            mnWav1 = single(mnWav1) * P.uV_per_bit;
        end
    end
    if ~isempty(vcFile)
        fclose(fid_bin);
    end
end %func
