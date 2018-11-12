%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function P = read_meta_file_(vcFile_meta)
    % Parse meta file, ask user if meta file doesn't exist
    P = [];
    try
        if exist(vcFile_meta, 'file') == 2
            S_meta = read_whisper_meta_(vcFile_meta);
            P = struct('sRateHz', S_meta.sRateHz, 'uV_per_bit', S_meta.scale, 'nChans', S_meta.nChans, 'vcDataType', S_meta.vcDataType);
            %'probe_file', [S_meta.vcProbe, '.prb'],
            P.Smeta = S_meta;
        else
            fprintf('%s is not found. Asking users to fill out the missing info\n', vcFile_meta);
            csAns = inputdlg_({...
            'sampling rate (Hz)', '# channels in file', ...
            'uV/bit', 'Header offset (bytes)', ...
            'Data Type (int16, uint16, single, double)', 'Neuropixels option (0 if N/A)'}, ...
            'Recording format', 1, {'30000', '385', '1','0','int16','0'});
            if isempty(csAns), return; end
            P = struct('sRateHz', str2double(csAns{1}), 'nChans', str2double(csAns{2}), ...
            'uV_per_bit', str2double(csAns{3}), 'header_offset', str2double(csAns{4}), ...
            'vcDataType', csAns{5}, 'imProbeOpt', str2double(csAns{6}));
            P.Smeta = P;
        end
    catch
        disperr_('read_meta_file_');
    end
end %func

%% local functions
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

%--------------------------------------------------------------------------
function S = imec3_imroTbl_(cSmeta)
    % Smeta has imroTbl

    vcDir_probe = 'C:\Dropbox (HHMI)\IMEC\SpikeGLX_Probe_Cal_Data\'; %this may crash. probe calibaration folder

    if isstruct(cSmeta), cSmeta = {cSmeta}; end %turn it into cell of struct
    % parse imroTbl
    cs_imroTbl = cellfun(@(S)S.imroTbl, cSmeta, 'UniformOutput', 0);
    cvn_imroTbl = cellfun(@(vc)textscan(vc, '%d', 'Delimiter', '( ),'), cs_imroTbl, 'UniformOutput', 0);
    cvn_imroTbl = cellfun(@(c)c{1}, cvn_imroTbl, 'UniformOutput', 0);
    S.viBank = cellfun(@(vn)vn(7), cvn_imroTbl);
    S.viRef = cellfun(@(vn)vn(8), cvn_imroTbl);
    S.vrGain_ap = single(cellfun(@(vn)vn(9), cvn_imroTbl));
    S.vrGain_lf = single(cellfun(@(vn)vn(10), cvn_imroTbl));
    S.nSites_bank = cvn_imroTbl{1}(4);

    Smeta1 = cSmeta{1};

    % correct gain
    nFiles = numel(S.viBank);
    nSites = numel(Smeta1.viSites);
    [mrScale_ap, mrScale_lf] = deal(ones(nSites, nFiles));
    S.vcProbeSN = sprintf('1%d%d', Smeta1.imProbeSN, Smeta1.imProbeOpt);
    % read gain correction
    vrGainCorr = ones(1, S.nSites_bank*4);
    if Smeta1.imProbeOpt ~= 2
        try
            vcFile_csv = sprintf('%s1%d%d\\Gain correction.csv', vcDir_probe, Smeta1.imProbeSN, Smeta1.imProbeOpt);
            try
                vrGainCorr = csvread(vcFile_csv, 1, 1);
            catch
                vrGainCorr = csvread(vcFile_csv, 0, 0);
            end
        catch
            ;
        end
    end

    % build scale
    for iFile = 1:nFiles
        vrGainCorr1 = vrGainCorr(Smeta1.viSites + double(S.nSites_bank*S.viBank(iFile)));
        mrScale_ap(:,iFile) = 1.2 * 1e6 / 2^10 / S.vrGain_ap(iFile) .* vrGainCorr1;
        mrScale_lf(:,iFile) = 1.2 * 1e6 / 2^10 / S.vrGain_lf(iFile) .* vrGainCorr1;
    end
    S.mrScale_ap = mrScale_ap;
    S.mrScale_lf = mrScale_lf;
end %func
