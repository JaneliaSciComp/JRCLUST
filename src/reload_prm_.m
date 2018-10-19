%--------------------------------------------------------------------------
function reload_prm_(hObject, event)
    % Edit prm file
    % 2016 07 06
    [S0, P] = get0_();
    S0.P = loadParam_(P.vcFile_prm);
    set(0, 'UserData', S0);
end %func
