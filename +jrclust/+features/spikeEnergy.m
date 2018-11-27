function en = spikeEnergy(spikeWindows)
    %SPIKEENERGY Summary of this function goes here
    %   Detailed explanation goes here
    en = shiftdim(std(spikeWindows, 1))';
end

