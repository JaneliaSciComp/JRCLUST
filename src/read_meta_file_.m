%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function P = read_meta_file_(vcFile_meta)
    %READ_META_FILE_ Read in a Whisper meta file (or ask for params for
    %some reason)

    P = [];
    try
        if exist(vcFile_meta, 'file') == 2
            S_meta = jrclust.utils.loadMetadata(vcFile_meta);
            P = struct('sRateHz', S_meta.sampleRate, ...
                       'uV_per_bit', S_meta.scale, ...
                       'nChans', S_meta.nChans, ...
                       'vcDataType', S_meta.dtype);
            P.Smeta = S_meta;
        else
            fprintf('%s is not found. Asking users to fill out the missing info\n', vcFile_meta);
            csAns = inputdlg_({...
                'sampling rate (Hz)', '# channels in file', ...
                'uV/bit', 'Header offset (bytes)', ...
                'Data Type (int16, uint16, single, double)', 'Neuropixels option (0 if N/A)'}, ...
                'Recording format', 1, {'30000', '385', '1','0','int16','0'});

            if isempty(csAns)
                return;
            end

            P = struct('sRateHz', str2double(csAns{1}), ...
                       'nChans', str2double(csAns{2}), ...
                       'uV_per_bit', str2double(csAns{3}), ...
                       'header_offset', str2double(csAns{4}), ...
                       'vcDataType', csAns{5}, 'imProbeOpt', str2double(csAns{6}));

            P.Smeta = P;
        end
    catch ME
        error('error reading meta file %s: %s', vcFile_meta, ME.message);
    end
end %func
