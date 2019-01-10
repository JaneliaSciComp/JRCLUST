function varargout = subsExt(filename, varargin)
    %SUBSEXT Substitute the extension of a file for a new extension
    %   [out1, out2, ..] = subsExt(filename, ext1, ext2, ...)
    [dirname, filename, ~] = fileparts(filename);
    if isempty(dirname)
        dirname = '.';
    end

    for i = 1:numel(varargin)
        newExt = varargin{i};
        varargout{i} = [dirname, filesep(), filename, newExt];
    end
end
