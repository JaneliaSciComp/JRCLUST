%--------------------------------------------------------------------------
function mr1 = fft_clean__(mr, thresh, nbins)
    % mr must be single

    if nargin<2, thresh = 6; end
    if nargin<3, nbins = 20; end
    nSkip_med = 4;
    nw = 3; %frequency neighbors to set to zero

    if thresh==0, thresh = []; end
    if isempty(thresh), mr1=mr; return ;end
    n = size(mr,1);
    n_pow2 = 2^nextpow2(n);
    vcClass = class_(mr);
    for iRetry = 1:2
        try
            mr1 = single(mr);
            vrMu = mean(mr1, 1);
            mr1 = bsxfun(@minus, mr1, vrMu);
            if n < n_pow2
                mr1 = fft(mr1, n_pow2);
            else
                mr1 = fft(mr1);
            end
            break;
        catch
            fprintf('GPU processing failed, retrying on CPU\n');
            mr = gather_(mr);
        end
    end %for

    % Find frequency outliers
    n1 = floor(n_pow2/2);
    viFreq = (1:n1)';
    % vrFft1 = abs(mean(bsxfun(@times, mr1(1+viFreq,:), viFreq), 2));
    vrFft1 = (mean(bsxfun(@times, abs(mr1(1+viFreq,:)), viFreq), 2));
    n2 = round(n1/nbins);
    for ibin=1:nbins
        vi1 = (n2*(ibin-1) : n2*ibin) + 1;
        if ibin==nbins, vi1(vi1>n1)=[]; end
        vrFft2 = vrFft1(vi1);
        vrFft2 = vrFft2 - median(vrFft2(1:nSkip_med:end)); %mad transform
        vrFft1(vi1) = vrFft2 / median(abs(vrFft2(1:nSkip_med:end)));
    end

    % broaden spectrum
    vl_noise = vrFft1>thresh;
    vi_noise = find(vl_noise);
    for i_nw=1:nw
        viA = vi_noise-i_nw;    viA(viA<1)=[];
        viB = vi_noise+i_nw;    viB(viB>n1)=[];
        vl_noise(viA)=1;
        vl_noise(viB)=1;
    end
    vi_noise = find(vl_noise);
    mr1(1+vi_noise,:) = 0;
    mr1(end-vi_noise+1,:) = 0;

    % inverse transform back to the time domain
    mr1 = real(ifft(mr1, n_pow2, 'symmetric')); %~30% faster than below
    if n < n_pow2, mr1 = mr1(1:n,:); end
    mr1 = bsxfun(@plus, mr1, vrMu); %add mean back
    mr1 = cast(mr1, vcClass); % cast back to the original type
end %func
