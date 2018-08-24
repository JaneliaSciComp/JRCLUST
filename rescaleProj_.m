%--------------------------------------------------------------------------
function rescaleProj_(vhPlot1, maxAmp, P)
    if nargin<3, P = get0_('P'); end
    for iPlot=1:numel(vhPlot1)
        hPlot1 = vhPlot1(iPlot);
        S_plot1 = get(hPlot1, 'UserData');
        update_plot2_proj_();
        S_plot1 = struct_delete_(S_plot1, 'hPoly'); %, 'hPlot_split'
        [vrX, vrY, viPlot, tr_dim] = amp2proj_(S_plot1.mrMin, S_plot1.mrMax, maxAmp, P.maxSite_show);
        S_plot1 = struct_add_(S_plot1, viPlot, vrX, vrY, maxAmp);
        set(hPlot1, 'XData', vrX, 'YData', vrY, 'UserData', S_plot1);
    end
end %func
