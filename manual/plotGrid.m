%--------------------------------------------------------------------------
function hPlot = plotGrid(lim, varargin)
    % lay out a lim(2) x lim(2) grid with tile boundaries at integer values

    vrX = floor((lim(1)*2:lim(2)*2 + 1)/2);
    vrY = repmat([lim(1), lim(2), lim(2), lim(1)], [1, ceil(numel(vrX)/4)]);
    vrY = vrY(1:numel(vrX));
    hPlot = plot([vrX(1:end-1), fliplr(vrY)], [vrY(1:end-1), fliplr(vrX)], varargin{:});
end % function
