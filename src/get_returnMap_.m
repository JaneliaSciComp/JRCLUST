%--------------------------------------------------------------------------
function [vrX, vrY] = get_returnMap_(iClu, P)
    vrTime1 = double(clu_time_(iClu)) / P.sRateHz;
    vrIsi1 = diff(vrTime1 * 1000); % in msec
    vrX = vrIsi1(1:end-1);
    vrY = vrIsi1(2:end);
    viShow = randperm(numel(vrX), min(P.nShow, numel(vrX)));
    vrX = vrX(viShow);
    vrY = vrY(viShow);
end
