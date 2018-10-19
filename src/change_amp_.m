%--------------------------------------------------------------------------
function [maxAmp, mrAmp_prev] = change_amp_(event, maxAmp, varargin)
    % varargin: plot object to rescale
    % Change amplitude scaling
    % change_amp_(event, maxAmp, varargin)
    % change_amp_(event) % directly set
    % if nargin<3, hPlot=[]; end
    factor = sqrt(2);
    if key_modifier_(event, 'shift'), factor = factor ^ 4; end
    mrAmp_prev = maxAmp;
    if strcmpi(event.Key, 'uparrow')
        maxAmp = maxAmp / factor;
    elseif strcmpi(event.Key, 'downarrow')
        maxAmp = maxAmp * factor;
    end
    for iPlot = 1:numel(varargin)
        try
            multiplot(varargin{iPlot}, maxAmp);
        catch
        end
    end
    % handle_fun_(@rescale_plot_, hPlot, maxAmp);
end %func
