%--------------------------------------------------------------------------
function S_mat = load_(vcFile, csVar, fVerbose)
    % return empty if the file doesn't exist
    % load_(vcFile)
    % load_(vcFile, vcVar)
    % load_(vcFile, csVar)
    % load_(vcFile, csVar, fVerbose)
    % load_(vcFile, [], fVerbose)

    S_mat = [];
    if nargin<2, csVar={}; end
    if ischar(csVar), csVar = {csVar}; end
    if nargin<3, fVerbose=1; end
    if exist(vcFile, 'file')
        try
            S_mat = load(vcFile);
        catch
            fprintf(2, 'Invalid .mat format: %s\n', vcFile);
        end
    else
        if fVerbose
            fprintf(2, 'File does not exist: %s\n', vcFile);
        end
    end
    if ~isempty(csVar)
        S_mat = get_(S_mat, csVar{:});
    end
end %func
