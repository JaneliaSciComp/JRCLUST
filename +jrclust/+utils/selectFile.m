function [filename, dirname] = selectFile(fileSpec, prompt, dirname, multiSelect)
    %SELECTFILE Select a file or files via a dialog
    args = {fileSpec, prompt, dirname};
    if multiSelect
        args = [args, {'MultiSelect', 'on'}];
    end

    [filename, dirname] = uigetfile(args{:});
    if ~(ischar(filename) || iscell(filename)) % cancel
        filename = jrclust.utils.ifEq(multiSelect, {''}, '');
        return;
    elseif ischar(filename) && multiSelect
        filename = {fullfile(dirname, filename)};
    elseif multiSelect
        filename = cellfun(@(f) fullfile(dirname, f), filename, 'UniformOutput', 0);
    end
end