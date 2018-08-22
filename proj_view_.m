%--------------------------------------------------------------------------
function proj_view_(hMenu)
    P = get0_('P');
    displayFeature = lower(get(hMenu, 'Label'));
    P.displayFeature = displayFeature;
    vhMenu = hMenu.Parent.Children;
    for iMenu=1:numel(vhMenu)
        vhMenu(iMenu).Checked = if_on_off_(vhMenu(iMenu).Label, displayFeature);
    end
    % auto-scale the view
    S0 = setUserData(P);
    button_CluWav_simulate_(S0.primarySelectedCluster, S0.secondarySelectedCluster, S0);
end %func
