%--------------------------------------------------------------------------
function figure_wait_(fWait, vhFig)
    % set all figures pointers to watch
    if nargin<2, vhFig = gcf; end
    if fWait
        set_(vhFig, 'Pointer', 'watch');
        drawnow;
    else
        set_(vhFig, 'Pointer', 'arrow');
    end
end %func
