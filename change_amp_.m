%--------------------------------------------------------------------------
function [maxAmp, prevMaxAmp] = change_amp_(event, maxAmp, varargin)
    % varargin: plot object to rescale
    % Change amplitude scaling
    % change_amp_(event, maxAmp, varargin)
    % change_amp_(event) % directly set

    factor = sqrt(2);

    if keyModifier(event, 'shift')
        factor = factor^4;
    end

    prevMaxAmp = maxAmp;

    if strcmpi(event.Key, 'uparrow')
        maxAmp = maxAmp/factor;
    elseif strcmpi(event.Key, 'downarrow')
        maxAmp = maxAmp*factor;
    end

    for iPlot = 1:numel(varargin)
        try
            multiplot(varargin{iPlot}, maxAmp);
        catch
        end
    end
end % function
