%--------------------------------------------------------------------------
function xylabel_(hAx, vcXLabel, vcYLabel, vcTitle)
    if isempty(hAx), hAx = gca; end
    xlabel(hAx, vcXLabel);
    ylabel(hAx, vcYLabel);
    if nargin>=4, title_(hAx, vcTitle); end
end %func
