%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_save_threshold_(hFig)
    % export S_fig.vnThresh and sets vcFile_thresh
    P = get0_('P');
    S_fig = get(hFig, 'UserData');
    vnThresh_site = S_fig.vnThresh_site;
    vcFile_thresh = strrep(P.prmFile, '.prm', '_thresh.mat');
    save(vcFile_thresh, 'vnThresh_site'); % also need to store filter values?
    P.vcFile_thresh = vcFile_thresh;
    set0_(P);
    edit_prm_file_(P, P.prmFile);

    msgbox_(sprintf('Saved to %s and updated %s (vcFile_thresh)', ...
    vcFile_thresh, P.prmFile));
end %func
