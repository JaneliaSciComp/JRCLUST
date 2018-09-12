%--------------------------------------------------------------------------
function [vi1, vrX1, vrY1, xlim1, ylim1] = imrect_plot_(hRect, hPlot1)
    % return points inside
    vrX = get(hPlot1, 'XData');
    vrY = get(hPlot1, 'YData');
    vrPos_rect = getPosition(hRect);
    xlim1 = vrPos_rect(1) + [0, vrPos_rect(3)];
    ylim1 = vrPos_rect(2) + [0, vrPos_rect(4)];
    vi1 = find(vrX>=xlim1(1) & vrX<=xlim1(2) & vrY>=ylim1(1) & vrY<=ylim1(2));
    if nargout>=2
        vrX1 = vrX(vi1);
        vrY1 = vrY(vi1);
    end
end %func
