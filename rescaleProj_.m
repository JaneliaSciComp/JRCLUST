%--------------------------------------------------------------------------
function rescaleProj_(hPlots, maxAmp, P)
    if nargin < 3
        P = get0_('P');
    end

    for iPlot = 1:numel(hPlots)
        hPlot = hPlots(iPlot);
        S_plot1 = get(hPlot, 'UserData');

        update_plot2_proj_();

        S_plot1 = struct_delete_(S_plot1, 'hPoly'); %, 'hPlot_split'
        
        if isempty(S_plot1)
            continue;
        end
        
        switch lower(P.displayFeature)
            case 'kilosort'
                % round to nearest hundred on either side of 0
                bounds = round(max(abs([S_plot1.mrMin(:) ; S_plot1.mrMax(:)]))/100, 0)*[-100, 100];
                maxPair = [];
                
            otherwise % vpp et al.
                bounds = [0 maxAmp];
                maxPair = P.maxSite_show;
        end

        [vrX, vrY, viPlot, ~] = amp2proj_(S_plot1.mrMin, S_plot1.mrMax, bounds, maxPair);

        S_plot1 = struct_add_(S_plot1, viPlot, vrX, vrY, maxAmp);
        set(hPlot, 'XData', vrX, 'YData', vrY, 'UserData', S_plot1);
    end
end %func
