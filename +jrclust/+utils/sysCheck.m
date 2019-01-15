function ok = sysCheck()
    %SYSCHECK Check for required toolboxes
    ok = 1;
    boxen = {'Signal', 'Image', 'Statistics'};
    reqs = cellfun(@(b) logical(license('test', [b '_Toolbox'])), boxen);

    if ~all(reqs)
        errmsg = ['You are missing the following toolboxes: ', strjoin(boxen(~reqs), ', ')]; 
        errordlg(errmsg, 'Missing system requirements');
        ok = 0;
    end
    
    if ~license('test', 'Distrib_Computing_Toolbox')
        wmsg = 'Could not find Parallel Toolbox (cannot use GPU)';
        warndlg(wmsg, 'Missing system requirements');
    end
end
