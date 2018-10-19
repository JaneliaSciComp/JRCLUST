%--------------------------------------------------------------------------
function close_hFig_traces_(hFig, event)
    try
        if ~ishandle(hFig), return; end
        if ~isvalid(hFig), return; end
        S_fig = get(hFig, 'UserData');
        fclose_(S_fig.fid_bin);
        try delete(hFig); catch; end %close one more time
    catch
        disperr_();
        close(hFig);
    end
end %func
