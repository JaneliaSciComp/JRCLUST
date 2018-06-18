%--------------------------------------------------------------------------
function sd2 = mr_std3_(mr, viRange1, viRange2)
    mr1 = mr(viRange1,:);
    mr2 = mr(viRange2,:);
    sd2 = sqrt(abs(mean(mr1.*mr2,2) - mean(mr1,2).*mean(mr2,2)));
end %func
