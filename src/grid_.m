%--------------------------------------------------------------------------
function grid_(vhAx, fGrid)
    if isempty(vhAx), vhAx = gca; end
    for iAx = 1:numel(vhAx)
        if ischar(fGrid)
            grid(vhAx(iAx), fGrid);
        else
            grid(vhAx(iAx), jrclust.utils.ifEq(fGrid, 'on', 'off'));
        end
    end
end %func
