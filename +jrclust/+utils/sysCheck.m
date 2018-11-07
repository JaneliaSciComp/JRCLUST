function sysCheck()
    % check for toolboxen
    toolBoxen = {'Signal', 'Image', 'Distrib_Computing', 'Statistics'};
    reqs = cellfun(@(b) logical(license('test', [b '_Toolbox'])), toolBoxen);
    if ~all(reqs)
        emsg = ['You are missing the following toolboxes: ', strjoin(toolBoxen(~reqs), ', ')]; 
        msgbox(emsg, 'Error', 'error');
        error(emsg);
    end
end
