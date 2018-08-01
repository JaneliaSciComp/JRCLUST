%--------------------------------------------------------------------------
function vhFig = get_fig_all_(csTag)
    vhFig = nan(size(csTag));
    for iFig=1:numel(csTag)
        try
            vhFig(iFig) = findobj('Tag', csTag{iFig});
        catch
            %         disperr_();
        end
    end %for
end %func
