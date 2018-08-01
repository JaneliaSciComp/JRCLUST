%--------------------------------------------------------------------------
function [mnWav2, vnWav2_mean] = filt_car_(mnWav2, P, mnWav1_pre, mnWav1_post, fTrim_pad)
    % Apply filter and CAR
    % @TODO: edge case
    if nargin<3, mnWav1_pre = []; end
    if nargin<4, mnWav1_post = []; end
    if nargin<5, fTrim_pad = 1; end
    n_pre = size(mnWav1_pre,1);
    n_post = size(mnWav1_post,1);
    if n_pre > 0 || n_post > 0
        mnWav2 = [mnWav1_pre; mnWav2; mnWav1_post];
    end
    P.vcFilter = get_filter_(P);
    switch lower(P.vcFilter)
        case 'user'
        %         vnFilter_user = -[5,0,-3,-4,-3,0,5]; % sgdiff acceleration
        vnFilter_user = single(get_set_(P, 'vnFilter_user', []));
        dialogAssert(~isempty(vnFilter_user), 'Set vnFilter_user to use vcFilter=''user''');
        for i=1:size(mnWav2,2)
            mnWav2(:,i) = conv(mnWav2(:,i), vnFilter_user, 'same');
        end
        case 'fir1'
        n5ms = round(P.sRateHz / 1000 * 5);
        vrFilter = single(fir1(n5ms, P.freqLim/P.sRateHz*2));
        for i=1:size(mnWav2,2)
            mnWav2(:,i) = conv(mnWav2(:,i), vrFilter, 'same');
        end
        case 'ndiff', mnWav2 = ndiff_(mnWav2, P.nDiff_filt);
        case 'fftdiff', mnWav2 = fftdiff_(mnWav2, P);
        %     case 'fftdiff', mnWav2 = fftdiff__(gather_(mnWav2), P.freqLim(2)/P.sRateHz/2);
        case {'sgdiff', 'sgfilt'}
        mnWav2 = sgfilt_(mnWav2, P.nDiff_filt);
        case 'bandpass'
        try
            mnWav2 = filtfilt_chain(single(mnWav2), P);
        catch
            fprintf('GPU filtering failed. Trying CPU filtering.\n');
            mnWav2 = filtfilt_chain(single(mnWav2), setfield(P, 'useGPU', 0));
        end
        mnWav2 = int16(mnWav2);
        case {'none', 'skip'} % no filter is applied
        ;
        case 'ndist'
        mnWav2 = ndist_filt_(mnWav2, get_set_(P, 'ndist_filt', 5));
        otherwise
        error('filt_car_: invalid filter option (vcFilter=''%s'')', P.vcFilter);
    end %switch

    % trim padding
    if (n_pre > 0 || n_post > 0) && fTrim_pad
        mnWav2 = mnWav2(n_pre+1:end-n_post,:);
    end

    %global subtraction before
    [mnWav2, vnWav2_mean] = wav_car_(mnWav2, P);
end %func
