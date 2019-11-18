function figPos = getDefaultFigPos(obj)
    figPos = containers.Map(); 
%     skipRD = ismember('FigRD',obj.hCfg.figList);
%     skipISI = ismember('FigISI',obj.hCfg.figList);
%     if ~skipRD && ~skipISI
%         figPos('FigTime') = [0.128 0.03 0.872 0.21];
%         figPos('FigHist') = [0.729 0.7 0.137 0.294];
%         figPos('FigSim') = [0.863 0.7 0.137 0.294];
%         figPos('FigCorr') = [0.596 0.7 0.137 0.294];
%         figPos('FigWav') = [0.128 0.24 0.47 0.755];
%         figPos('FigProj') = [0.596 0.24 0.404 0.46];
%         figPos('FigMap') = [0 0.61 0.13 0.385];
%         figPos('FigPos') = [0 0.03 0.13 0.58];
%     elseif ~skipISI
%         figPos('FigCorr') = [0.584 0.7 0.209 0.295];
%         figPos('FigSim') = [0.791 0.7 0.209 0.295];
%         figPos('FigWav') = [0.128 0.24 0.458 0.755];
%         figPos('FigHist') = [0.85 0.355 0.15 0.344];
%         figPos('FigProj') = [0.584 0.24 0.268 0.46];
%         figPos('FigISI') = [0.85 0.03 0.149 0.325];
%         figPos('FigTime') = [0.128 0.03 0.724 0.21];
%         figPos('FigPos') = [0 0.03 0.13 0.58];
%         figPos('FigMap') = [0 0.61 0.13 0.385];
%     else
%         figPos('FigTime') = [0.128 0.03 0.72 0.21];
%         figPos('FigRD') = [0.701 0.7 0.147 0.295];
%         figPos('FigWav') = [0.128 0.24 0.429 0.755];
%         figPos('FigSim') = [0.847 0.7 0.153 0.295];
%         figPos('FigProj') = [0.555 0.24 0.293 0.46];
%         figPos('FigHist') = [0.847 0.355 0.153 0.344];
%         figPos('FigISI') = [0.847 0.03 0.153 0.325];
%         figPos('FigCorr') = [0.555 0.7 0.147 0.295];
%         figPos('FigPos') = [0 0.03 0.13 0.58];
%         figPos('FigMap') = [0 0.61 0.13 0.385];   
%     end
    figPos('FigPos') = [0 0 .15 .5];
    figPos('FigMap') = [0 .5 .15 .5];
    figPos('FigWav') = [.15 .2 .35 .8];
    figPos('FigProj') = [.5 .2 .35 .5];
    figPos('FigSim') = [.5 .7 .35 .3];
    figPos('FigRD') = [.85 0 .15 .25];    
    skipRD = ~ismember('FigRD', obj.hCfg.figList); % skip the rho-delta plot
    if skipRD % expand figTime et al. to take its place
        figPos('FigTime') = [.15 0 .85 .2];
        figPos('FigCorr') = [.85 .2 .15 .27];
        figPos('FigISI') = [.85 .47 .15 .26];
        figPos('FigHist') = [.85 .73 .15 .27];
    else
        figPos('FigTime') = [.15 0 .7 .2];
        figPos('FigCorr') = [.85 .25 .15 .25];
        figPos('FigISI') = [.85 .5 .15 .25];
        figPos('FigHist') = [.85 .75 .15 .25];
    end
end