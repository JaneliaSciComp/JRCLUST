%--------------------------------------------------------------------------
function [mrFet1, mrFet2, mrFet3, trWav2_spk] = trWav2fet_(tnWav1_spk, P, nSites_spk, spikeSecondarySites)
    % [mrFet1, mrFet2, mrFet3, trWav_spk2] = trWav2fet_(tnWav_spk1, P)
    % mrFet = trWav2fet_(tnWav_spk1, P)
    if nargin<3, nSites_spk = []; end
    if nargin<4, spikeSecondarySites = []; end

    [mrFet1, mrFet2, mrFet3] = deal(single([]));
    trWav2_spk = single(permute(tnWav1_spk, [1,3,2]));
    trWav2_spk = spkwav_car_(trWav2_spk, P, nSites_spk, spikeSecondarySites);
    % if get_set_(P, 'fMeanSubt_fet', 1), trWav2_spk = meanSubt_(trWav2_spk); end % 12/16/17 JJJ experimental

    switch lower(P.vcFet) %{'xcor', 'amp', 'slope', 'pca', 'energy', 'vpp', 'diff248', 'spacetime'}
        case {'spacetime', 'cov', 'cov2'}
        [mrFet1, mrFet2] = trWav2fet_cov_(trWav2_spk, P);
        case 'cov_prev'
        nDelay = 3;
        gtrWav1 = meanSubt_(trWav2_spk);
        mr1 = zscore_(gtrWav1(:,:,1));
        mr2 = zscore_(gtrWav1([ones(1,nDelay),1:end-nDelay],:,1));
        mrFet1 = mean(gtrWav1 .* repmat(mr1, [1,1,size(gtrWav1,3)]), 1);
        mrFet2 = mean(gtrWav1 .* repmat(mr2, [1,1,size(gtrWav1,3)]), 1);
        mrFet1 = shiftdim(mrFet1,1)';
        mrFet2 = shiftdim(mrFet2,1)';

        case {'vpp', 'vppsqrt'}
        mrFet1 = shiftdim(max(trWav2_spk) - min(trWav2_spk))';
        if strcmpi(P.vcFet, 'vppsqrt'), mrFet1 = sqrt(mrFet1); end

        case {'amp', 'vmin'}
        mrFet1 = shiftdim(abs(min(trWav2_spk)))';

        case {'vminmax', 'minmax'}
        mrFet1 = shiftdim(abs(min(trWav2_spk)))';
        mrFet2 = shiftdim(abs(max(trWav2_spk)))';

        case 'energy'
        mrFet1 = shiftdim(std(trWav2_spk,1))';

        case 'energy2'
        nDelay = 3;
        mrFet1 = shiftdim(std(trWav2_spk,1))';
        trcov_ = @(a,b)shiftdim(sqrt(abs(mean(a.*b) - mean(a).*mean(b))));
        mrFet2 = trcov_(trWav2_spk(1:end-nDelay,:,:), trWav2_spk(nDelay+1:end,:,:))';
        %mrFet1 = shiftdim(std(trWav_spk1,1))';

        case {'pca', 'gpca', 'fpca'}
        %  Compute PrinVec, 2D, max channel only
        if strcmpi(P.vcFet, 'fpca')
            trWav2_spk0 = trWav2_spk;
            trWav2_spk = fft(trWav2_spk);
            trWav2_spk = abs(trWav2_spk(1:end/2,:,:));
        end
        if strcmpi(P.vcFet, 'pca')
            mrPv = tnWav2pv_(trWav2_spk, P);
        else
            %             trWav_spk1 = spkwav_car_(trWav_spk1, viSites_ref);
            mrPv_global = get0_('mrPv_global');
            if isempty(mrPv_global)
                [mrPv_global, vrD_global] = tnWav2pv_(trWav2_spk, P);
                [mrPv_global, vrD_global] = gather_(mrPv_global, vrD_global);
                setUserData(mrPv_global, vrD_global);
            end
            mrPv = mrPv_global;
        end
        [mrFet1, mrFet2, mrFet3] = project_interp_(trWav2_spk, mrPv, P);
    end
    if nargout==1
        switch P.nPcPerChan
            case 2
            mrFet1 = cat(1, mrFet1, mrFet2);
            case 3
            mrFet1 = cat(1, mrFet1, mrFet2, mrFet3);
        end %switch
    end
end %func
