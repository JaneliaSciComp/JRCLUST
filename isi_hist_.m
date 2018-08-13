%--------------------------------------------------------------------------
function vnHist = isi_hist_(iClu1, vrX)
    P = get0_('P');
    vrTime1 = double(clu_time_(iClu1)) / P.sampleRateHz;
    vnHist = hist(diff(vrTime1)*1000, vrX);
    vnHist(end)=0;
    vnHist = vnHist ./ sum(vnHist);
end
