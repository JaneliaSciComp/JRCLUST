%--------------------------------------------------------------------------
function batch_(vcFile_batch, vcCommand)
    % batch process parameter files (.batch) file
    % batch_(myfile.batch, vcCommand): contains a list of .prm files
    % batch_(vcFile_batch, vcFile_template): batch contains .bin files

    if nargin<2, vcCommand=[]; end
    if isempty(vcCommand), vcCommand = 'spikesort'; end

    if matchFileExt_(vcCommand, '.prm')
        vcFile_template = vcCommand;
        vcCommand = 'spikesort';
    else
        vcFile_template = default_prm_path_();
    end

    csFiles_prm = load_batch_(vcFile_batch);
    vlFile_err = true(size(csFiles_prm));
    csErr_file = cell(size(csFiles_prm));

    for iFile = 1:numel(csFiles_prm)
        vcFile_prm_ = csFiles_prm{iFile};
        if matchFileExt_(vcFile_prm_, {'.bin', '.dat'});
            vcFile_prm_ = makeprm_(vcFile_prm_, vcFile_template, 0);
        end
        try
            jrc('clear');
            fprintf('Processing file %d/%d: %s\n', iFile, numel(csFiles_prm), vcFile_prm_);
            jrc(vcCommand, vcFile_prm_);
            vlFile_err(iFile) = 0;
        catch
            csErr_file{iFile} = lasterr();
            fprintf(2, 'Error in file %d/%d: %s\n', iFile, numel(csFiles_prm), vcFile_prm_);
            fprintf(2, '\t%s\n\n', csErr_file{iFile});
        end
    end %for

    % summary
    if any(vlFile_err) %error found
        fprintf(2,'Errors found in the .prm files below:\n');
        viFile_err = find(vlFile_err);
        for iErr = 1:numel(viFile_err)
            vcFile_err_ = csFiles_prm{viFile_err(iErr)};
            fprintf(2,'\tError %d/%d: %s\n', iErr, numel(viFile_err), vcFile_err_);
            fprintf(2,'\t\t%s\n\n', csErr_file{iFile});
        end
    else
        fprintf('All %d files processed successfully\n', numel(csFiles_prm));
    end
    edit(vcFile_batch);
end %func
