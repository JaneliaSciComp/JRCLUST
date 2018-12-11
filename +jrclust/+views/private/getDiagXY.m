function [xVals, yVals] = getDiagXY(lim)
    %GETDIAGXY Get x-y values for outlines of diagonal entries
    %   lim: [start, end] or [start, end, offset]

    yVals_ = floor((lim(1)*2:lim(2)*2 + 1)/2);
    xVals_ = [yVals_(2:end), lim(2)];

    xVals = [xVals_(1:end - 1), fliplr(yVals_)];
    yVals = [yVals_(1:end - 1), fliplr(xVals_)];

    if numel(lim) >= 3 % offset
        xVals = xVals + lim(3);
        yVals = yVals + lim(3);
    end
end
