%--------------------------------------------------------------------------
function hPlot = plotTable_(lim, varargin)

    XData = floor((lim(1)*2:lim(2)*2+1)/2);
    YData = repmat([lim(1), lim(2), lim(2), lim(1)], [1, ceil(numel(XData)/4)]);
    YData = YData(1:numel(XData));
    hPlot = plot([XData(1:end-1), fliplr(YData)], [YData(1:end-1), fliplr(XData)], varargin{:});
end %func
