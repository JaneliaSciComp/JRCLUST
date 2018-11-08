function ok = sysCheck()
    %SYSCHECK Check for required toolboxes
    ok = true;
    boxen = {'Signal', 'Image', 'Distrib_Computing', 'Statistics'};
    reqs = cellfun(@(b) logical(license('test', [b '_Toolbox'])), boxen);

    if ~all(reqs)
        emsg = ['You are missing the following toolboxes: ', strjoin(boxen(~reqs), ', ')]; 
        errordlg(emsg, 'Missing system requirements');
        ok = false;
    end
end
