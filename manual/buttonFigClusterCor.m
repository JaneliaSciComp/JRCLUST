%--------------------------------------------------------------------------
function buttonFigClusterCor(xyPos, vcButton)
    S0 = get(0, 'UserData');
    xyPos = round(xyPos);

    if xyPos(1) < 1
        xyPos(1) = 1;
    end
    if xyPos(2) < 1
        xyPos(2) = 1;
    end

    switch lower(vcButton)
        case 'normal' %left click
            S0.primarySelectedCluster = xyPos(1);

            if diff(xyPos) == 0
                S0.secondarySelectedCluster = [];
            else
                S0.secondarySelectedCluster = xyPos(2);
            end

            S0 = button_CluWav_simulate_(S0.primarySelectedCluster, S0.secondarySelectedCluster, S0);
            S0 = keyPressFcn_cell_(getCachedFig('FigWav'), {'z'}, S0); %zoom
    end %switch
end % function
