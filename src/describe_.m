%--------------------------------------------------------------------------
function csDesc = describe_(vcFile_prm)
    %describe_()
    %describe_(vcFile_prm)
    %describe_(S0)

    % describe _jrc.mat file
    if nargin==0
        S0 = get(0, 'UserData');
    elseif isstruct(vcFile_prm)
        S0 = vcFile_prm;
    else
        S0 = load(strrep(vcFile_prm, '.prm', '_jrc.mat'));
    end
    P = S0.P;

    nSites = numel(P.viSite2Chan);
    tDur = double(max(S0.viTime_spk) - min(S0.viTime_spk)) / P.sRateHz;
    nSpk = numel(S0.viTime_spk);
    nSitesPerEvent = P.maxSite*2+1;


    csDesc = {};
    csDesc{end+1} = sprintf('Recording file');
    csDesc{end+1} = sprintf('    Recording file          %s', P.vcFile);
    csDesc{end+1} = sprintf('    Probe file              %s', P.probe_file);
    csDesc{end+1} = sprintf('    Recording Duration      %0.1fs ', tDur);
    csDesc{end+1} = sprintf('    #Sites                  %d', nSites);
    csDesc{end+1} = sprintf('Events');
    csDesc{end+1} = sprintf('    #Spikes                 %d', nSpk);
    csDesc{end+1} = sprintf('    Feature                 %s', P.vcFet);
    csDesc{end+1} = sprintf('    #Sites/event            %d', nSitesPerEvent);
    if ~isempty(get_(S0, 'nLoads'))
        csDesc{end+1} = sprintf('    #Loads                  %d', S0.nLoads);
    end

    if isfield(S0, 'S_clu')
        S_clu = S0.S_clu;
        csDesc{end+1} = sprintf('Cluster');
        csDesc{end+1} = sprintf('    #Clusters               %d', S_clu.nClu);
        csDesc{end+1} = sprintf('    #Unique events          %d', sum(S_clu.viClu>0));
        csDesc{end+1} = sprintf('    min. spk/clu            %d', P.min_count);
        if isfield(S_clu, 't_runtime')
            csDesc{end+1} = sprintf('    Cluster run-time        %0.1fs', S_clu.t_runtime);
        end
    end
    try
        runtime_total = S0.runtime_detect + S0.runtime_sort;
        csDesc{end+1} = sprintf('Runtime (s)');
        csDesc{end+1} = sprintf('    Detect + feature        %0.1fs', S0.runtime_detect);
        csDesc{end+1} = sprintf('    Cluster                 %0.1fs', S0.runtime_sort);
        csDesc{end+1} = sprintf('    Total                   %0.1fs', runtime_total);
        csDesc{end+1} = sprintf('    Runtime speed           x%0.1f realtime', tDur / runtime_total);
    catch
        ;
    end

    if nargout==0
        cellfun(@(x)disp(x), csDesc);
    end

end %func
