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
    if ~isfield(P, 'displayFeature'), P.displayFeature = 'vpp'; end
    [vrFet1, viSpk1] = getFet_clu_(iClu, iSite, S0);
    vrTime1 = double(S0.spikeTimes(viSpk1)) / P.sampleRateHz;

    % label
    switch lower(P.displayFeature)
        case {'vpp', 'vmin'} %voltage feature
        vcYlabel = sprintf('Site %d (\\mu%s)', iSite, P.displayFeature);
        otherwise %other feature options
        vcYlabel = sprintf('Site %d (%s)', iSite, P.displayFeature);
    end

end %func
