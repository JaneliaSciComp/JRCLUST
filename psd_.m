%--------------------------------------------------------------------------
function [mrPower, vrFreq] = psd_(mr, Fs, nSkip)
    % power spectrum
    if nargin<3, nSkip = 1; end
    if nSkip>1
        mr = mr(1:nSkip:end,:);
        %     Fs = Fs / nSkip;
    end
    n = size(mr,1);
    n1 = round(n/2);
    mrPower = fft(meanSubtract(single(mr)));
    mrPower = pow2db_(abs(mrPower(2:n1+1, :))) / n;

    if nargout>=2, vrFreq = Fs*(1:n1)'/n; end
end % function
