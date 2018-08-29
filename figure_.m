%--------------------------------------------------------------------------
% error tolerent figure selection
function hFig = figure_(hFig)
    if nargin<1, hFig = figure(); return; end
    if isempty(hFig), return; end
    try
        if gcf() ~= hFig
            figure(hFig);
        end
    catch; end
end % function
