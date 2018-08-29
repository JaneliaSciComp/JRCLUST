%--------------------------------------------------------------------------
function [mrPc1, mrPc2, mrPv1, mrPv2] = pca_pc_spk_(viSpk1, viSites1, mrPv1, mrPv2)
    % varargin: mrPv
    % varargout: mrPc
    % project
    %[mrPc1, mrPc2, mrPv1, mrPv2] = pca_pc_spk_(viSpk1, viSites1)
    %[mrPc1, mrPc2] = pca_pc_spk_(viSpk1, viSites1, mrPv1, mrPv2)
    nSites1 = numel(viSites1);
    spikeWaveforms1 = permute(getSpikeWaveformsSites(viSpk1, viSites1, [], 0), [1,3,2]);
    if nargin < 3
        [mrPv1, mrPv2] = pca_pv_spk_(viSpk1, viSites1, spikeWaveforms1);
    end

    dimm1 = size(spikeWaveforms1); %nT x nSpk x nChan
    [mrPc1, mrPc2] = deal(zeros(dimm1(2), nSites1, 'single'));

    try
        for iSite1=1:nSites1
            mrWav_spk1 = meanSubtract(single(spikeWaveforms1(:,:,iSite1)));
            mrPc1(:,iSite1) = (mrPv1(:,iSite1)' * mrWav_spk1)';
            mrPc2(:,iSite1) = (mrPv2(:,iSite1)' * mrWav_spk1)';
        end %for
    catch
        disperr_();
    end

    mrPc1 = (mrPc1') / dimm1(1);
    mrPc2 = (mrPc2') / dimm1(1);
end % function
