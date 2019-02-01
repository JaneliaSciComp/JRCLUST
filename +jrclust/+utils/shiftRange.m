function [shiftsBack, shiftsForward] = shiftRange(nSamples, nShifts, shiftFactors)
    %SHIFTRANGE
    if ~isempty(nShifts)
        [shiftsBack, shiftsForward] = deal(cell(nShifts*2 + 1, 1));
        shiftFactors = -nShifts:nShifts;
    else
        [shiftsBack, shiftsForward] = deal(cell(numel(shiftFactors), 1));
    end

    sampleRange = 1:nSamples;
    for iFactor = 1:numel(shiftFactors)
        iShift = shiftFactors(iFactor);

        back = -round(iShift/2);
        fwd  = iShift + back;

        rangeBack = sampleRange + back;
        rangeForward = sampleRange + fwd;
        mask = (rangeBack >= 1 & rangeBack <= nSamples) & (rangeForward >=1 & rangeForward <= nSamples);

        shiftsBack{iFactor} = rangeBack(mask);
        shiftsForward{iFactor} = rangeForward(mask);
    end
end
