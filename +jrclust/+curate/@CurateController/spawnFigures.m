function spawnFigures(obj)
    %SPAWNFIGURES Create the standard cadre of figures
    obj.hFigs = containers.Map();   
    for f=1:length(obj.hCfg.figList)
        figTag = obj.hCfg.figList{f};
        figToolbar = 0;
        figMenubar = 0;
        switch figTag
            case 'FigPos'
                figTitle = 'Unit position';
                figToolbar = 1;
            case 'FigMap'
                figTitle = 'Probe map';
                figToolbar = 1;
            case 'FigWav'
                figTitle = 'Averaged waveform';
                figMenubar = 1;
            case 'FigTime'
                figTitle = 'Feature vs. time';
            case 'FigProj'
                figTitle = 'Feature projection';
            case 'FigSim'
                figTitle = 'Template-based similarity score';                                    
            case 'FigHist'
                figTitle = 'ISI histogram';                
            case 'FigISI'
                figTitle = 'Return map';                
            case 'FigCorr'
                figTitle = 'Time correlation';                
            case 'FigRD'
                if isa(obj.hClust,'jrclust.sort.TemplateClustering')
                    warning('Skipping spawning of rho-delta plot because density-peak clustering was not used.');
                    continue
                end
                figTitle = 'Unit rho-delta';                
        end
        obj.hFigs(figTag) = jrclust.views.Figure(figTag,obj.hCfg.figPos{f},sprintf('%s: %s',figTitle,obj.hCfg.sessionName),figToolbar,figMenubar);
        drawnow;
    end
end
