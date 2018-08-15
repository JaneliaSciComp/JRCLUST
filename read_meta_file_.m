%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function P = read_meta_file_(vcFile_meta)
    % Parse meta file, ask user if meta file doesn't exist
    P = [];
    try
        if exist(vcFile_meta, 'file') == 2
            S_meta = read_whisper_meta_(vcFile_meta);
            P = struct('sampleRateHz', S_meta.sampleRateHz, 'uV_per_bit', S_meta.scale, 'nChans', S_meta.nChans, 'dataType', S_meta.dataType);
            %'probeFile', [S_meta.vcProbe, '.prb'],
            P.Smeta = S_meta;
        else
            fprintf('%s is not found. Asking users to fill out the missing info\n', vcFile_meta);
            csAns = inputdlg_({...
            'sampling rate (Hz)', '# channels in file', ...
            'uV/bit', 'Header offset (bytes)', ...
            'Data Type (int16, uint16, single, double)', 'Neuropixels option (0 if N/A)'}, ...
            'Recording format', 1, {'30000', '385', '1','0','int16','0'});
            if isempty(csAns), return; end
            P = struct('sampleRateHz', str2double(csAns{1}), 'nChans', str2double(csAns{2}), ...
            'uV_per_bit', str2double(csAns{3}), 'headerOffset', str2double(csAns{4}), ...
            'dataType', csAns{5}, 'imProbeOpt', str2double(csAns{6}));
            P.Smeta = P;
        end
    catch
        disperr_('read_meta_file_');
    end
end %func
