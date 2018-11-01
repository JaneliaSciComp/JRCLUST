%--------------------------------------------------------------------------
function [vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, iClu, S0)
    % just specify iSite to obtain background info
    % 2016 07 07 JJJ
    % return feature correspojnding to a site and cluster
    % requiring subsampled info: cvrVpp_site and cmrFet_site. store in S0

    if nargin < 2, iClu = []; end
    if nargin<3, S0 = get(0, 'UserData'); end
    % S_clu = S0.S_clu;
    P = S0.P;
    if ~isfield(P, 'vcFet_show'), P.vcFet_show = 'vpp'; end
    [vrFet1, viSpk1] = getFet_clu_(iClu, iSite, S0);
    vrTime1 = double(S0.viTime_spk(viSpk1)) / P.sRateHz;

    % label
    switch lower(P.vcFet_show)
        case {'vpp', 'vmin'} %voltage feature
            vcYlabel = sprintf('Site %d (\\mu%s)', iSite, P.vcFet_show);

        otherwise %other feature options
            vcYlabel = sprintf('Site %d (%s)', iSite, P.vcFet_show);
    end

end %func
