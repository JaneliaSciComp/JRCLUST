%--------------------------------------------------------------------------
% 17/11/10: Removed recursive saving
function S0 = load0_(vcFile_mat)
    % Load a mat file structure and set to 0 structure
    % S0 = load0_(vcFile_mat)
    % S0 = load0_(P)
    % only set the S0 if it's found
    if isstruct(vcFile_mat)
        P = vcFile_mat;
        vcFile_mat = strrep(P.vcFile_prm, '.prm', '_jrc.mat');
    end

    if ~exist(vcFile_mat, 'file')
        S0 = [];
        fprintf('File %s does not exist\n', vcFile_mat);
        return;
    end

    try
        fprintf('loading %s...\n', vcFile_mat); t1=tic;
        S0 = load(vcFile_mat);
        set(0, 'UserData', S0);
        fprintf('\ttook %0.1fs\n', toc(t1));
    catch
        S0 = [];
        disperr_();
    end
    if isfield(S0, 'S0'), S0 = rmfield(S0, 'S0'); end % Remove recursive saving
end %func
