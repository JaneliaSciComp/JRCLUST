%--------------------------------------------------------------------------
% 12/20/17 JJJ: Export to LFP file
function import_lfp_(P)
    % % Merge LFP file for IMEC3 probe
    % try
    %     if ~isempty(strfind(lower(vcFile), '.imec.ap.bin'))
    %         func_ap2lf = @(x)strrep(lower(x), '.imec.ap.bin', '.imec.lf.bin');
    %         vcFile_lf = func_ap2lf(vcFile);
    %         multiFilenames_lf = cellfun(@(x)func_ap2lf(x), multiFilenames1, 'UniformOutput', 0);
    %         merge_binfile_(vcFile_lf, multiFilenames_lf);
    %     end
    % catch
    %     disp('Merge LFP file error for IMEC3.');
    % end
    P.lfpFile = strrep(P.paramFile, '.prm', '.lfp.jrc');
    t1 = tic;
    if ~isfield(P, 'multiFilenames') || isempty(P.multiFilenames)
        % single file
        if is_new_imec_(P.vcFile) % don't do anything, just set the file name
            P.lfpFile = strrep(P.vcFile, '.imec.ap.bin', '.imec.lf.bin');
            P.nSkip_lfp = 12;
            P.sampleRateHz_lfp = 2500;
        else
            bin_file_copy_(P.vcFile, P.lfpFile, P);
        end
    else % craete a merged output file
        csFiles_bin = filter_files_(P.multiFilenames);
        fid_lfp = fopen(P.lfpFile, 'w');
        % multiple files merging
        P_ = P;
        for iFile = 1:numel(csFiles_bin)
            vcFile_ = csFiles_bin{iFile};
            if is_new_imec_(vcFile_)
                vcFile_ = strrep(vcFile_, '.imec.ap.bin', '.imec.lf.bin');
                P_.nSkip_lfp = 1;
            end
            bin_file_copy_(vcFile_, fid_lfp, P_);
        end
        fclose(fid_lfp);
    end
    % update the lfp file name in the parameter file
    updateParamFile(P, P.paramFile);
    fprintf('\tLFP file (lfpFile) updated: %s\n\ttook %0.1fs\n', P.lfpFile, toc(t1));
end % function
