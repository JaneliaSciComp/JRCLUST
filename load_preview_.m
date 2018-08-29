%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function [mnWav_raw, S_preview] = load_preview_(P)
    % Load the subsampled dataset
    % Useful for inspecting threshold and so on. filter and
    % S_preview: which file and where it came from
    if ischar(P), P = loadParams(P); end

    nLoads_max_preview = getOr(P, 'nLoads_max_preview', 30);
    sec_per_load_preview = getOr(P, 'sec_per_load_preview', 1);

    % determine files to load
    if ~isfield(P, 'multiFilenames') || isempty(P.multiFilenames)
        csFile_bin = {P.vcFile};
    else
        csFile_bin = filter_files_(P.multiFilenames);
    end
    csFile_bin = subsample_(csFile_bin, nLoads_max_preview);

    % load files
    nLoads_per_file = floor(nLoads_max_preview / numel(csFile_bin));
    nSamples_per_load = round(sec_per_load_preview * P.sampleRateHz);

    % file loading loop
    [mnWav_raw, cviLim_load, csFile_load] = deal({});
    % [mnWav_raw, mnWav_filt] = deal({});
    P.useGPU = 0;
    for iFile = 1:numel(csFile_bin)
        try
            vcFile_bin_ = csFile_bin{iFile};
            [fid_bin_, nBytes_bin_, P.headerOffset] = fopenInfo(vcFile_bin_, 'r');
            setUserData(P);
            if isempty(fid_bin_)
                fprintf(2, '.bin file does not exist: %s\n', vcFile_bin_);
                continue;
            end
            fprintf('File %d/%d: %s\n\t', iFile, numel(csFile_bin), vcFile_bin_);
            nSamples_bin_ = floor(nBytes_bin_ / bytesPerSample_(P.dataType) / P.nChans);
            if nSamples_bin_ < nSamples_per_load % load the whole thing
                nLoads_per_file_ = 1;
                nSamples_per_load_ = nSamples_bin_;
            else
                nLoads_per_file_ = min(nLoads_per_file, floor(nSamples_bin_ / nSamples_per_load));
                nSamples_per_load_ = nSamples_per_load;
            end
            [cvi_lim_bin, viRange_bin] = sample_skip_([1, nSamples_per_load_], nSamples_bin_, nLoads_per_file_);
            for iLoad_ = 1:nLoads_per_file_
                fprintf('.');
                ilim_bin_ = cvi_lim_bin{iLoad_};
                fseek_(fid_bin_, ilim_bin_(1), P);
                mnWav_raw{end+1} = load_file_(fid_bin_, diff(ilim_bin_) + 1, P);
                cviLim_load{end+1} = ilim_bin_;
                csFile_load{end+1} = vcFile_bin_;
            end
            fprintf('\n');
        catch
            fprintf(2, 'Loading error: %s\n', vcFile_bin_);
        end
        fclose_(fid_bin_, 0);
    end
    nLoads = numel(mnWav_raw);
    mnWav_raw = cell2mat(mnWav_raw');
    % if nargout>=2, mnWav_raw = cell2mat(mnWav_raw'); end
    if nargout>=2
        S_preview = makeStruct_(nLoads_per_file, nLoads_max_preview, ...
        sec_per_load_preview, nSamples_per_load, nLoads, csFile_load, cviLim_load);
    end
end % function
