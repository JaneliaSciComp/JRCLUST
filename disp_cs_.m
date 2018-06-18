%--------------------------------------------------------------------------
function disp_cs_(cs)
    % display cell string
    cellfun(@(s)fprintf('%s\n',s), cs);
end %func
