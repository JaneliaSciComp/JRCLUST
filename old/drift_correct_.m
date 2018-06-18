%--------------------------------------------------------------------------
function [mrFet, vrY_pre, vrY_post] = drift_correct_(mrFet, P) % drift correction
    % Manually correct drift using {[t1,t2,d1,d2], ...} format in cvrDepth_drift
    % mrFet last column contains the y coordinate
    fPlot_drift = 1;
    if isempty(getfield_(P, 'cvrDepth_drift')), return; end
    S0 = get(0, 'UserData');
    fprintf('Correcting drift\n\t'); t1=tic;
    vrY_pre = mrFet(end,:);
    vrY_post = vrY_pre;
    for iDrift = 1:numel(P.cvrDepth_drift)
        tlim_dlim1 = P.cvrDepth_drift{iDrift};
        tlim1 = tlim_dlim1(1:2) * P.sRateHz;
        dlim1 = tlim_dlim1(3:4);
        viSpk1 = find(S0.viTime_spk >= tlim1(1) & S0.viTime_spk < tlim1(2));
        if isempty(viSpk1), continue; end
        viTime_spk1 = S0.viTime_spk(viSpk1);
        if diff(dlim1) == 0
            vrY_post(viSpk1) = vrY_pre(viSpk1) + dlim1(1); % use single depth correction factor
        else
            vrY_post(viSpk1) = vrY_pre(viSpk1) + interp1(tlim1, dlim1, viTime_spk1, 'linear', 'extrap');
        end
        fprintf('.');
    end
    mrFet(end,:) = vrY_post;
    fprintf('\n\tDrift correction took %0.1fs\n', toc(t1));

    if fPlot_drift
        [viTime_spk, vrAmp_spk] = get0_('viTime_spk', 'vrAmp_spk');
        viSpk = find(vrAmp_spk < median(vrAmp_spk)); %pick more negative
        %     viSpk1 = viSpk1(1:2:end); %plot every other
        figure; hold on; set(gcf,'Color','w');
        ax(1)=subplot(121);
        plot(viTime_spk(viSpk), vrY_pre(viSpk), 'go', 'MarkerSize', 2); title('before correction'); grid on;
        ax(2)=subplot(122);
        plot(viTime_spk(viSpk), vrY_post(viSpk), 'bo', 'MarkerSize', 2); title('after correction'); grid on;
        linkaxes(ax,'xy');
    end
end %func
