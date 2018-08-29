%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function gui_(vcArg1, vcFile_prm_)
    % JRCLUST GUI interface

    if ~isempty(vcArg1)
        vcFile_prm = vcArg1;
    elseif ~isempty(vcFile_prm_)
        vcFile_prm = vcFile_prm_;
    else
        vcFile_prm = '';
    end
    S_gui.vcFile_prm = vcFile_prm;
    setUserData(S_gui);
    jrc_gui(vcFile_prm);

end % function
