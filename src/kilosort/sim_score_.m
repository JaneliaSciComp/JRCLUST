%--------------------------------------------------------------------------
function S_clu = sim_score_(S_clu)
    %SIM_SCORE_ update the pairwise similarity scores of clusters
    %

    S0 = get(0, 'UserData');
    P = S0.P;

    if get_set_(P, 'fImportKsort', 0)
        rez = S0.rez;

        nClu = S_clu.nClu;

        viTemplate_spk = S_clu.viTemplate_spk;
        viClu = S_clu.viClu;

        % update the unique template indices for each cluster
        S_clu.cviTemplate_clu = arrayfun(@(iClu) unique(viTemplate_spk(viClu == iClu)), 1:S_clu.nClu, 'UniformOutput', 0);

        mrSim_clu = zeros(nClu);
        for iClu = 1:nClu
            viTemp_clu = S_clu.cviTemplate_clu{iClu}; % unique template indices for spikes in this cluster

            % compute cluster sim score, Phy style
            sims = max(rez.simScore(viTemp_clu, :), [], 1);

            for jClu = iClu:nClu
                viTemp_clu2 = S_clu.cviTemplate_clu{jClu};
                mrSim_clu(iClu, jClu) = max(sims(viTemp_clu2));
                mrSim_clu(jClu, iClu) = mrSim_clu(iClu, jClu);
            end
        end

        S_clu.mrSim_clu = mrSim_clu;
    end
end
