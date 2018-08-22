%--------------------------------------------------------------------------
function S0 = button_CluWav_simulate_(primarySelectedCluster, secondarySelectedCluster, S0, initial_load) % TW
    if nargin < 4
        initial_load = 0;
    end

    if nargin<3
        S0 = get(0, 'UserData');
    end

    if nargin<2
        secondarySelectedCluster = [];
    end

    if primarySelectedCluster == secondarySelectedCluster
        secondarySelectedCluster = [];
    end

    hFig = gcf;
    figure_wait_(1, hFig);

    S0 = update_cursor_(S0, primarySelectedCluster, 0);
    S0 = update_cursor_(S0, secondarySelectedCluster, 1);
    if initial_load==0
        S0 = keyPressFcn_cell_(getCachedFig('FigWav'), {'t','c','i','v','e','f'}, S0); %'z' to recenter 'j', %TW
    else
        S0 = keyPressFcn_cell_(getCachedFig('FigWav'), {'j','t','c','i','v','e','f'}, S0); %'z' to recenter 'j', %TW
    end
    set(0, 'UserData', S0);

    plot_raster_(S0); %psth
    figure_wait_(0, hFig);
end
