%--------------------------------------------------------------------------
function [vrX, vrY] = plotDiag__(lim)
    % lim: [start, end] or [start, end, offset]

    vrY0 = floor((lim(1)*2:lim(2)*2+1)/2);
    % vrY0 = lim(1):lim(end);
    vrX0 = [vrY0(2:end), lim(2)];
    vrX = [vrX0(1:end-1), fliplr(vrY0)];
    vrY = [vrY0(1:end-1), fliplr(vrX0)];
    if numel(lim)>=3
        vrX = vrX + lim(3);
        vrY = vrY + lim(3);
    end
end %func
