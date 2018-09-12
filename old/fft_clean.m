function mr1 = fft_clean(mr, thresh, nbins)
% mr must be single

if nargin<2, thresh = 6; end
if nargin<3, nbins = 20; end
nSkip_med = 4;
nw = 3; %frequency neighbors to set to zero
fDebug = 0;

if thresh==0, thresh = []; end
if isempty(thresh), mr1=mr; return ;end

% mr = single(mr);
vrMu = mean(mr, 1);
mr1 = bsxfun(@minus, mr, vrMu);
n = size(mr1,1);
n_pow2 = 2^nextpow2(n);
if n < n_pow2
    mr1 = fft(mr1, n_pow2);
else
    mr1 = fft(mr1);
end
n1 = floor(n_pow2/2);
viFreq = (1:n1)';
% vrFft1 = abs(mean(bsxfun(@times, mr1(1+viFreq,:), viFreq), 2));
vrFft1 = (mean(bsxfun(@times, abs(mr1(1+viFreq,:)), viFreq), 2));
if fDebug || nargout==0 
    vrFft0 = mean(abs(mr1(1+viFreq,:)), 2); vrFreq = (1:numel(vrFft1))/numel(vrFft1)*20000/2;  
    figure; subplot(221); plot(vrFreq, 2*pow2db(vrFft0),'k.','MarkerSize',8); xlabel('Freq (Hz)'); ylabel('Power (dB)'); grid on; xlim([0 500]);     
end

n2 = round(n1/nbins); 
for ibin=1:nbins
    vi1 = (n2*(ibin-1) : n2*ibin) + 1;
    if ibin==nbins, vi1(vi1>n1)=[]; end
    vrFft2 = vrFft1(vi1);
%     vrFft2 = detrend(vrFft2);
    vrFft2 = vrFft2 - median(vrFft2(1:nSkip_med:end)); %mad transform
    vrFft1(vi1) = vrFft2 / median(abs(vrFft2(1:nSkip_med:end)));
%     vrFft1(vi1) = zscore((vrFft1(vi1)));    
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
if fDebug || nargout==0 
    subplot(222); plot(vrFreq, vrFft1,'k.','MarkerSize',8); xlabel('Freq (Hz)'); ylabel('z-score (detrended)'); 
    grid on; xlim([0 500]); ylim([0 50]); set(gca,'YScale', 'linear');
    hold on; plot(get(gca,'XLim'),thresh*[1 1], 'r-','MarkerSize',8);
    hold on; plot(vrFreq(vi_noise), vrFft1(vi_noise), 'r.','MarkerSize',8);
    % hold on; plot(vrFreq(vi_noise), vrPowDb(vi_noise), 'r.');
end

mr1(1+vi_noise,:) = 0;
mr1(end-vi_noise+1,:) = 0;
mr1 = real(ifft(mr1, n_pow2, 'symmetric')); %~30% faster than below
% mr1 = real(ifft(mr1));
if n < n_pow2, mr1 = mr1(1:n,:); end
mr1 = bsxfun(@plus, mr1, vrMu); %add mean back

if nargout==0 || fDebug
    figure; set(gcf,'Color','w');
    subplot(121); plot(vrFreq, vrFft1,'k.','MarkerSize',8); xlabel('Freq (Hz)');
    axis([0 10000 -10 40]);
    n = min(10000, size(mr,1));
    subplot(122); plot(mr(1:n,1)); hold on; plot(mr1(1:n,1));
%     plot(mr1(1:n,1)-mr(1:n,1));    
    
%     return;
%     sRateHz = 25000;
%     vrFreq = (1:size(mr1,1))/size(mr1,1)*sRateHz;
%     figure; plot(viFreq, vrFft0, '.');
    subplot(1,2,2); plot(viFreq, vrFft1, '.');
    hold on; plot(viFreq(vi_noise), vrFft1(vi_noise), 'o');
    grid on; ylim([-5 20]);
    
    plot(get(gca,'XLim'), thresh*[1,1]);
    axis([0 10000 -10 40]);
end