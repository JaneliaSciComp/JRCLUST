%--------------------------------------------------------------------------
function proj_view_(hMenu)
    P = get0_('P');
    vcFet_show = lower(get(hMenu, 'Label'));
    P.vcFet_show = vcFet_show;
    vhMenu = hMenu.Parent.Children;
    for iMenu=1:numel(vhMenu)
        vhMenu(iMenu).Checked = if_on_off_(vhMenu(iMenu).Label, vcFet_show);
    end
    % auto-scale the view
    S0 = set0_(P);
    button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
end %func
