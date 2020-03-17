function figPos = getCurrentFigPos()
    % utility that reports the position of the current figure windows
    figs = get(0,'Children');
    figPos = containers.Map();
    for f=1:length(figs)
        figs(f).Units = 'normalized';
        figPos(figs(f).Tag) = round(figs(f).OuterPosition,3);
        fprintf('figPos(''%s'') = [%g %g %g %g]\n',figs(f).Tag,round(figs(f).OuterPosition,3));
    end
end