%--------------------------------------------------------------------------
function hPlot = plotDiag_(lim, varargin)
    [vrX, vrY] = plotDiag__(lim);
    % vrY = floor((lim(1)*2:lim(2)*2+1)/2);
    % vrX = [vrY(2:end), lim(end)];
    % hPlot = plot([vrX(1:end-1), fliplr(vrY)], [vrY(1:end-1), fliplr(vrX)], varargin{:});
    hPlot = plot(vrX, vrY, varargin{:});
end %func
