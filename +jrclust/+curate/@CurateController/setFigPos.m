function setFigPos(figPos)
    % utility for manually setting the position of the current figure windows
    figs = get(0,'Children');
    keys = figPos.keys;
    tags = {figs.Tag};
    for f=1:length(keys)
        idx = strcmp(tags,keys{f});
        if any(idx)
            figs(find(idx)).Units = 'normalized';
            figs(find(idx)).OuterPosition = figPos(keys{f});
        else
            warning('%s not a current figure tag.',keys{f});
        end
    end
end