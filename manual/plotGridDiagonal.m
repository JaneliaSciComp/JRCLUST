%--------------------------------------------------------------------------
function hPlot = plotGridDiagonal(lim, varargin)
    [vrX, vrY] = getGridDiagonal(lim);
    hPlot = plot(vrX, vrY, varargin{:});
end % function
