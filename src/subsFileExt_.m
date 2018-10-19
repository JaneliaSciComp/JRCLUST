%--------------------------------------------------------------------------
% 8/2/17 JJJ: added '.' if dir is empty
% 7/31/17 JJJ: Substitute file extension
function varargout = subsFileExt_(vcFile, varargin)
    % Substitute the extension part of the file
    % [out1, out2, ..] = subsFileExt_(filename, ext1, ext2, ...)

    [vcDir_, vcFile_, ~] = fileparts(vcFile);
    if isempty(vcDir_), vcDir_ = '.'; end
    for i=1:numel(varargin)
        vcExt_ = varargin{i};
        varargout{i} = [vcDir_, filesep(), vcFile_, vcExt_];
    end
end %func
