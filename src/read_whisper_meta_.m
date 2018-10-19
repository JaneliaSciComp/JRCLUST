%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function S = read_whisper_meta_(vcFname)
    % Import SpikeGLX meta file format

    S = [];
    viRef_imec3 = [37 76 113 152 189 228 265 304 341 380];

    % Read file
    if nargin < 1
        [FileName,PathName,FilterIndex] = uigetfile();
        vcFname = fullfile(PathName, FileName);
        if ~FilterIndex
            return;
        end
    end

    try
        %Read Meta
        S = text2struct_(vcFname);
        S.vcDataType = 'int16'; %whisper standard
        S.ADC_bits = 16;
        S.vcProbe = '';

        %convert new fields to old fields
        if isfield(S, 'niSampRate')
            % SpikeGLX
            S.nChans = S.nSavedChans;
            S.sRateHz = S.niSampRate;
            S.rangeMax = S.niAiRangeMax;
            S.rangeMin = S.niAiRangeMin;
            S.auxGain = S.niMNGain;
            try
                S.outputFile = S.fileName;
                S.sha1 = S.fileSHA1;
                S.vcProbe = 'imec2';
            catch
                S.outputFile = '';
                S.sha1 = [];
            end
        elseif isfield(S, 'imSampRate')
            % IMECIII
            S.nChans = S.nSavedChans;
            S.sRateHz = S.imSampRate;
            S.rangeMax = S.imAiRangeMax;
            S.rangeMin = S.imAiRangeMin;
            S.ADC_bits = 10; %10 bit adc but 16 bit saved
            vnIMRO = textscan(S.imroTbl, '%d', 'Delimiter', '( ),');
            vnIMRO = vnIMRO{1};
            S.auxGain = double(vnIMRO(9)); %hard code for now;
            S.auxGain_lfp = double(vnIMRO(10)); %hard code for now;
            S.vcProbe = sprintf('imec3_opt%d', vnIMRO(3));
            S.nSites = vnIMRO(4);
            S.viSites = setdiff(1:S.nSites, viRef_imec3); %sites saved
            try
                S.S_imec3 = imec3_imroTbl_(S);
            catch
                S.S_imec3 = [];
            end
        elseif isfield(S, 'sample_rate') %nick steinmetz
            S.nChans = S.n_channels_dat;
            S.sRateHz = S.sample_rate;
        end

        %number of bits of ADC [was 16 in Chongxi original]
        try
            S.scale = ((S.rangeMax-S.rangeMin)/(2^S.ADC_bits))/S.auxGain * 1e6; %uVolts
        catch
            S.scale = 1;
        end
        S.uV_per_bit = S.scale;
        if isfield(S, 'auxGain_lfp')
            S.uV_per_bit_lfp = S.scale * S.auxGain / S.auxGain_lfp;
        end
    catch
        disp(lasterr);
    end
end %func
