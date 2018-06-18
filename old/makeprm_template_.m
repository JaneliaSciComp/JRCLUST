%--------------------------------------------------------------------------
function vcFile_prm = makeprm_template_(vcFile_bin, vcFile_template, vcFile_prb)
    % output prm file from a template file
    % vcFile_prm is [vcFile_bin, vcFile_prb, '.prm'], after removing .bin and .prb extensions
    csLines_prm = {};
    csLines_prm{end+1} = sprintf('vcFile = ''%s'';', vcFile_bin);
    csLines_prm{end+1} = sprintf('template_file = ''%s'';', vcFile_template);
    if ~isempty(vcFile_prb)
        csLines_prm{end+1} = sprintf('probe_file = ''%s'';', vcFile_prb);
    else
        S_prm = file2struct_(vcFile_template);
        vcFile_prb = S_prm.probe_file;
    end

    % update from meta file if exists
    S_meta = read_meta_file_(subsFileExt_(vcFile_bin, '.meta'));
    if ~isempty(S_meta)
        csLines_prm{end+1} = sprintf('nChans = %d;', S_meta.nChans);
        csLines_prm{end+1} = sprintf('sRateHz = %d;', S_meta.sRateHz);
        csLines_prm{end+1} = sprintf('uV_per_bit = %f;', S_meta.uV_per_bit);
    end

    % name prm file
    [~,vcPostfix,~] = fileparts(vcFile_prb);
    vcFile_prm = subsFileExt_(vcFile_bin, ['_', vcPostfix, '.prm']);

    cellstr2file_(vcFile_prm, csLines_prm);
end
