%--------------------------------------------------------------------------
function export_rate_()
    S_clu = get0_('S_clu');
    mrRate_clu = clu_rate_(S_clu);
    csMsg = assignWorkspace_(mrRate_clu);
    fprintf(csMsg);
    jrclust.utils.qMsgBox(csMsg, 1);
end %func
