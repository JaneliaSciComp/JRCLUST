%--------------------------------------------------------------------------
function [viSpk, vrSpk, viSite] = spikeMerge_(cviSpk, cvrSpk, P)
    % provide spike index (cviSpk) and amplitudes (cvrSPk) per sites

    nSites = numel(cviSpk);
    viSpk = jrclust.utils.neCell2mat(cviSpk);      vrSpk = jrclust.utils.neCell2mat(cvrSpk);
    viSite = jrclust.utils.neCell2mat(cellfun(@(vi,i)repmat(i,size(vi)), cviSpk, num2cell((1:nSites)'), 'UniformOutput', false));
    [viSpk, viSrt] = sort(viSpk);   vrSpk = vrSpk(viSrt);   viSite = viSite(viSrt);
    viSite = int32(viSite);
    viSpk = int32(viSpk);

    [cviSpkA, cvrSpkA, cviSiteA] = deal(cell(nSites,1));

    try
        parfor iSite = 1:nSites %parfor speedup: 2x %parfor
            try
                [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P);
            catch
                disperr_();
            end
        end
    catch
        for iSite = 1:nSites
            try
                [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P);
            catch
                disperr_();
            end
        end
    end

    % merge parfor output and sort
    viSpk = jrclust.utils.neCell2mat(cviSpkA);
    vrSpk = jrclust.utils.neCell2mat(cvrSpkA);
    viSite = jrclust.utils.neCell2mat(cviSiteA);
    [viSpk, viSrt] = sort(viSpk); %sort by time
    vrSpk = jrclust.utils.tryGather(vrSpk(viSrt));
    viSite = viSite(viSrt);
end %func
