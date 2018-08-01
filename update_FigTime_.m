%--------------------------------------------------------------------------
% function [vrX1_spk, vrY1_spk, vrVpp_spk, mrVpp1, trWav1] = spikePos_(viSpk1, viSite1, P)
% % Get spike position from a spike index viSpk1
% % Determine Vpp
% % viTime1 = randomSelect(viTime1, P.nShow_proj); %select fewer
% % [trWav1, mrVpp1] = mr2tr_spk_(mrWav, viTime1, viSite1, P);
%
% fCentroid = 0;
%
% if fCentroid
%     P.fInterp_mean = 0;
% %     [mrWav_mean1, imax_site1, trWav_int1] = interp_align_mean_(trWav1, P);
%     mrXY_spk = centroid_pca_(trWav1, P.mrSiteXY(viSite1, :));
% end
% vrX1_site = P.mrSiteXY(viSite1, 1);
% vrY1_site = P.mrSiteXY(viSite1, 2);
% mrVpp1_sq = mrVpp1.^2;
% vrVpp1_sq_sum = sum(mrVpp1_sq);
% if ~fCentroid
%     vrX1_spk = sum(bsxfun(@times, mrVpp1_sq, vrX1_site)) ./ vrVpp1_sq_sum;
%     vrY1_spk = sum(bsxfun(@times, mrVpp1_sq, vrY1_site)) ./ vrVpp1_sq_sum;
% else
%     vrX1_spk = mrXY_spk(:,1);
%     vrY1_spk = mrXY_spk(:,1);
% end
% vrVpp_spk = sqrt(vrVpp1_sq_sum);
% end %func


%--------------------------------------------------------------------------
function update_FigTime_()
    % display features in a new site

    [hFig, S_fig] = get_fig_cache_('FigTime');
    S0 = get(0, 'UserData');
    P = S0.P;
    if ~isVisible_(S_fig.hAx), return ;end
    % P.vcFet_show = S_fig.csFet{S_fig.iFet};
    setUserData(P);
    [vrFet0, vrTime0, vcYlabel] = getFet_site_(S_fig.iSite, [], S0);
    if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
    toggleVisible_(S_fig.hPlot0, S_fig.fPlot0);
    update_plot_(S_fig.hPlot0, vrTime0, vrFet0);
    set(S_fig.hPlot1, 'YData', getFet_site_(S_fig.iSite, S0.iCluCopy, S0));
    % set(S_fig.hPlot0, 'XData', vrTime0, 'YData', vrFet0);
    if ~isempty(S0.iCluPaste)
        set(S_fig.hPlot2, 'YData', getFet_site_(S_fig.iSite, S0.iCluPaste, S0));
    else
        hide_plot_(S_fig.hPlot2);
        %     set(S_fig.hPlot2, 'XData', nan, 'YData', nan);
    end
    % switch lower(P.vcFet_show)
    %     case 'vpp'
    ylim_(S_fig.hAx, [0, 1] * S_fig.maxAmp);
    imrect_set_(S_fig.hRect, [], [0, 1] * S_fig.maxAmp);
    %     otherwise
    %         ylim_(S_fig.hAx, [0, 1] * P.maxAmp);
    %         imrect_set_(S_fig.hRect, [], [0, 1] * P.maxAmp);
    % end
    grid(S_fig.hAx, 'on');
    ylabel(S_fig.hAx, vcYlabel);
end %func
