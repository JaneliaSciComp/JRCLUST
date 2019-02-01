function exportToWorkspace(vals, verbose)
    %EXPORTTOWORKSPACE Export a values in a struct to local workspace
    fieldNames = fieldnames(vals);
    varsAssigned = cell(numel(fieldNames), 1);

    for i = 1:numel(fieldNames)
        fn = fieldNames{i};
        if ~isempty(vals.(fn))
            assignin('base', fn, vals.(fn)); % TODO: check for overwrite
            varsAssigned{i} = fn;
        end
    end

    if verbose
        varsAssigned = varsAssigned(cellfun(@(va) ~isempty(va), varsAssigned));
        msg = 'The following have been assigned to workspace:';
        msg = sprintf('%s\n\t* %s', msg, strjoin(varsAssigned, '\n\t* '));
        jrclust.utils.qMsgBox(msg);
    end
end
