function S = loadMetadata(metafile)
    %LOADMETADATA convert SpikeGLX metadata to struct
    %   TODO: build into Config workflow to cut down on specified parameters

    % get absolute path of metafile
    metafile_ = jrclust.utils.absPath(metafile);
    if isempty(metafile_)
        error('could not find meta file %s', metafile);
    end

    try
        S = jrclust.utils.metaToStruct(metafile_);
    catch ME
        error('could not read meta file %s: %s', metafile, ME.message);
    end

    S.adcBits = 16;
    S.probe = '';
    S.isImec = 0;

    %convert new fields to old fields
    if isfield(S, 'niSampRate') % SpikeGLX
        S.nChans = S.nSavedChans;
        S.sampleRate = S.niSampRate;
        S.rangeMax = S.niAiRangeMax;
        S.rangeMin = S.niAiRangeMin;
        S.gain = S.niMNGain;
        try
            S.outputFile = S.fileName;
            S.sha1 = S.fileSHA1;
            S.probe = 'imec2';
        catch
            S.outputFile = '';
            S.sha1 = [];
        end
    elseif isfield(S, 'imSampRate') % IMEC probe
        S.nChans = S.nSavedChans;
        S.sampleRate = S.imSampRate;
        S.rangeMax = S.imAiRangeMax;
        S.rangeMin = S.imAiRangeMin;
        
        S.probeOpt = [];
        % Determine probe type: 3A (0), 3B (1), or NP2.0 (2)
        if isfield(S,'imProbeOpt')
            probeType = '3A';
            S.probeOpt = S.imProbeOpt;
        elseif isfield(S,'imDatPrb_type')
            if S.imDatPrb_type == 0
                probeType = '3B'; 
            elseif S.imDatPrb_type == 21 || S.imDatPrb_type == 24                
                probeType = 'NP2';
            else
                probeType = 'unknown';
            end
        else
            probeType = 'unknown';
        end
        
        if strcmp(probeType,'3A') || strcmp(probeType, '3B')
            % 3A or 3B data; both have 10 bit adc, gain specified in imro
            S.adcBits = 10; % 10 bit adc but 16 bit saved
            % read data from ~imroTbl
            imroTbl = strsplit(S.imroTbl(2:end-1), ')(');
            % parse first channel entry
            imroTblChan = cellfun(@str2double, strsplit(imroTbl{2}, ' '));
            S.gain = imroTblChan(4);
            S.gainLFP = imroTblChan(5);
        elseif strcmp(probeType,'NP2')
            % NP 2.0 -- headstage has two docks          
            S.adcBits = 14; % 14 bit adc but 16 bit saved
            S.gain = 80; % constant gain
        end

        S.isImec = 1;
    end

    %number of bits of ADC [was 16 in Chongxi original]
    try
        S.scale = ((S.rangeMax - S.rangeMin)/(2^S.adcBits))/S.gain * 1e6; %uVolts
    catch
        S.scale = 1;
    end

    S.bitScaling = S.scale;
    if isfield(S, 'gainLFP')
        S.bitScalingLFP = S.scale * S.gain / S.gainLFP;
    end
end
