function validateFigs(obj)
    %% validates figure list and positions, and fills in with defaults if no positions specified
    if ~isempty(obj.figPos) 
        if numel(obj.figList)~=numel(obj.figPos)
            error('Error parsing params: size of figList does not match size of figPos in param file!');
        end
    else %% use default figure positions
        obj.figPos = defaultFigPos(obj.figList);
    end

    %% add back in required figs if they are missing
    requiredFigs = {'FigWav'};
    if ~all(ismember(requiredFigs,obj.figList))
        warnMsg = sprintf('"%s" will be used despite absence from user-defined figure list.\n',requiredFigs{:}); 
        warndlg(warnMsg,'Missing required figure in "FigList" param');
        for i = 1:length(requiredFigs)
            if ~ismember(requiredFigs{i},obj.figList)
                obj.figList = [obj.figList requiredFigs(i)];
                obj.figPos = [obj.figPos defaultFigPos(requiredFigs{i})];
            end
        end
    end
end

function figPos = defaultFigPos(figList)
    figPos = cell(1, numel(figList));
    hasFigRD = ismember('FigRD', figList);

    for f=1:length(figList)
       switch figList{f}
           case 'FigCorr'
               if hasFigRD
                   figPos{f} = [.85 .25 .15 .25];
               else
                   figPos{f} = [.85 .2 .15 .27];
               end
               
           case 'FigHist'
               if hasFigRD
                   figPos{f} = [.85 .75 .15 .25];
               else
                   figPos{f} = [.85 .73 .15 .27];
               end

           case 'FigISI'
               if hasFigRD
                   figPos{f} = [.85 .5 .15 .25];
               else
                   figPos{f} = [.85 .47 .15 .26];
               end

           case 'FigMap'
               figPos{f} = [0 .5 .15 .5];

           case 'FigPos'
               figPos{f} = [0 0 .15 .5];

           case 'FigProj'
               figPos{f} = [.5 .2 .35 .5];

           case 'FigRD'
               figPos{f} = [.85 0 .15 .25];

           case 'FigSim'
               figPos{f} = [.5 .7 .35 .3];

           case 'FigTime'
               if hasFigRD
                   figPos{f} = [.15 0 .7 .2];
               else
                   figPos{f} = [.15 0 .85 .2];
               end

           case 'FigWav'
               figPos{f} = [.15 .2 .35 .8];
       end
    end
end
