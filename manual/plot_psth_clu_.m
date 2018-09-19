%--------------------------------------------------------------------------
function plot_psth_clu_(viTime_clu, vrTime_trial, P, hAx, vcColor)
    if nargin<4, hAx=gca; end
    if nargin<5, vcColor = 'k'; end

    tbin = P.tbin_psth;
    nbin = round(tbin * P.sRateHz);
    nlim = round(P.tlim_psth/tbin);
    viTime_Trial = round(vrTime_trial / tbin);

    vlTime1=zeros(0);
    vlTime1(ceil(double(viTime_clu)/nbin))=1;
    mr1 = vr2mr2_(double(vlTime1), viTime_Trial, nlim);
    vnRate = mean(mr1,2) / tbin;
    vrTimePlot = (nlim(1):nlim(end))*tbin + tbin/2;
    bar(hAx, vrTimePlot, vnRate, 1, 'EdgeColor', 'none', 'FaceColor', vcColor);
    vrXTick = P.tlim_psth(1):(P.xtick_psth):P.tlim_psth(2);
    set(hAx, 'XTick', vrXTick, 'XTickLabel', []);
    grid(hAx, 'on');
    hold(hAx, 'on');
    plot(hAx, [0 0], get(hAx,'YLim'), 'r-');
    ylabel(hAx, 'Rate (Hz)');
    xlim_(hAx, P.tlim_psth);
end
