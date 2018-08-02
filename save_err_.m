%--------------------------------------------------------------------------
% 17/12/5 JJJ: created
function save_err_(hErr, vcMsg)
    % save error object
    if nargin<2, vcMsg = ''; end

    vcDatestr = datestr(now, 'yyyy_mmdd_HHMMSS');
    S_err = struct('message', hErr.message(), 'error', vcMsg, ...
    'P', get0_('P'), 'time', vcDatestr, 'MException', hErr);

    vcVar_err = sprintf('err_%s', vcDatestr);
    eval(sprintf('%s = S_err;', vcVar_err));

    % save to a file
    vcFile_err = fullfile(jrcpath_(), 'jrc_log_error.mat');
    if ~fileExists(vcFile_err)
        save(vcFile_err, vcVar_err, '-v7.3');
    else
        try
            save(vcFile_err, vcVar_err ,'-append');
        catch
            save(vcFile_err, vcVar_err, '-v7.3'); % file might get too full, start over
        end
    end
    fprintf('Error log updated: %s\n', vcFile_err);
end %func
