%--------------------------------------------------------------------------
function export_rate_()
    S_clu = get0_('S_clu');
    mrRate_clu = clu_rate_(S_clu);
    csMsg = assignWorkspace_(mrRate_clu);
    fprintf(csMsg);
    msgbox_(csMsg, 1);
end %func
