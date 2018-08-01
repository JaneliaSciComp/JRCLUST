%--------------------------------------------------------------------------
function update_plot2_proj_(vrX, vrY)
    if nargin==0, vrX=nan; vrY=nan; end
    [hFig, S_fig] = getCachedFig('FigProj');
    % erase polygon
    if nargin==0
        try
            update_plot_(S_fig.hPlot2, vrX, vrY);
            delete(findobj(get(S_fig.hAx, 'Child'), 'Type', 'hggroup'));
        catch
            ;
        end
    end
end
