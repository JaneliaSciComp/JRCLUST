%--------------------------------------------------------------------------
function [xvals, yvals, viPlot, tr_dim] = amp2proj_(yvals, xvals, bounds, maxPair, P)
    if nargin < 4
        maxPair = [];
    end
    if nargin < 5
        P = get0_('P');
    end

    xvals = linmap_(xvals', bounds, [0,1], 1);
    yvals = linmap_(yvals', bounds, [0,1], 1);

    [nSpikes, nChans] = size(yvals);
    if isempty(maxPair)
        maxPair = nChans;
    end

    % spike features translated into site-site boxes
    [transX, transY] = deal(nan([nSpikes, nChans, nChans], 'single'));

    for jChan = 1:nChans
        jsiteY = yvals(:, jChan);
        ymask = jsiteY > 0  & jsiteY < 1; % get points away from the boundaries
        for iChan = 1:nChans
            if abs(iChan - jChan) > maxPair
                continue;
            end
            
            % vpp only:
            % min on site j vs. min on site i above the diagonal
            if strcmpi(P.vcFet_show, 'vpp') && jChan > iChan
                isiteX = yvals(:, iChan);
            else % diagonal and below: min vs. max
                isiteX = xvals(:, iChan);
            end

            xymask = (isiteX > 0 & isiteX < 1) & ymask;

            transX(xymask, jChan, iChan) = isiteX(xymask) + iChan - 1;
            transY(xymask, jChan, iChan) = jsiteY(xymask) + jChan - 1;
        end
    end

    % plot projection
    viPlot = find(~isnan(transX) & ~isnan(transY));
    xvals = transX(viPlot); xvals = xvals(:);
    yvals = transY(viPlot); yvals = yvals(:);
    tr_dim = size(transX);
end %func

%% local functions
function vals = linmap_(vals, prevLim, newLim, saturate)
    if nargin < 4
        saturate = 0;
    end

    if numel(prevLim) == 1
        prevLim = abs(prevLim)*[-1, 1];
    end
    if numel(newLim) == 1
        newLim = abs(newLim)*[-1, 1];
    end

    if saturate
        vals(vals > prevLim(2)) = prevLim(2);
        vals(vals < prevLim(1)) = prevLim(1);
    end

    if prevLim(1) == prevLim(2) % % ignore newLim and just rescale 
        vals = vals / prevLim(1);
    else
        vals = interp1(prevLim, newLim, vals, 'linear', 'extrap');
    end
end %func
