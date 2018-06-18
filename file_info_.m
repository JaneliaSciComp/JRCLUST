%--------------------------------------------------------------------------
% 7/31/17 JJJ: documentation and testing
function P = file_info_(vcFile)
    % Returns empty if file not found or multiple ones found

    S_dir = dir(vcFile);
    if numel(S_dir)==1
        P = struct('vcDate_file', S_dir.date, 'nBytes_file', S_dir.bytes);
    else
        P = [];
    end
end %func
