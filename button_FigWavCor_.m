%--------------------------------------------------------------------------
function button_FigWavCor_(xyPos, vcButton)
    S0 = get(0, 'UserData');
    xyPos = round(xyPos);
    switch lower(vcButton)
        case 'normal' %left click
        S0.iCluCopy = xyPos(1);
        if diff(xyPos) == 0
            S0.iCluPaste = [];
        else
            S0.iCluPaste = xyPos(2);
        end
        S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
        S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
    end %switch
end %func
