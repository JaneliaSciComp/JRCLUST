%--------------------------------------------------------------------------
function [filteredTraces, vnWav2_mean] = filt_car_(rawTraces, P, prePadding, postPadding, fTrim_pad)
    % Apply filter and CAR
    % @TODO: edge case

    if nargin < 3
        prePadding = [];
    end
    if nargin < 4
        postPadding = [];
    end
    if nargin < 5
        fTrim_pad = 1;
    end

    n_pre = size(prePadding, 1);
    n_post = size(postPadding, 1);

    filteredTraces = [prePadding; rawTraces; postPadding];
    P.vcFilter = get_filter_(P);

    switch lower(P.vcFilter)
        case 'user'
            vnFilter_user = single(getOr(P, 'vnFilter_user', []));
            dialogAssert(~isempty(vnFilter_user), 'Set vnFilter_user to use vcFilter=''user''');
            for i = 1:size(filteredTraces, 2)
                filteredTraces(:, i) = conv(filteredTraces(:, i), vnFilter_user, 'same');
            end
        case 'fir1'
            n5ms = round(P.sampleRateHz / 1000 * 5);
            vrFilter = single(fir1(n5ms, P.freqLim/P.sampleRateHz*2));
            for i = 1:size(filteredTraces,2)
                filteredTraces(:,i) = conv(filteredTraces(:,i), vrFilter, 'same');
            end
        case 'ndiff'
            filteredTraces = ndiff_(filteredTraces, P.nDiff_filt);
        case 'fftdiff'
            filteredTraces = fftdiff_(filteredTraces, P);
        case {'sgdiff', 'sgfilt'}
            filteredTraces = sgfilt_(filteredTraces, P.nDiff_filt);
        case 'bandpass'
            try
                filteredTraces = filtfilt_chain(single(filteredTraces), P);
            catch
                fprintf('GPU filtering failed. Trying CPU filtering.\n');
                filteredTraces = filtfilt_chain(single(gather(filteredTraces)), setfield(P, 'useGPU', 0));
            end
            filteredTraces = int16(filteredTraces);
        case 'ndist'
            filteredTraces = ndist_filt_(filteredTraces, getOr(P, 'ndist_filt', 5));
        case {'none', 'skip'}
            ; % nothing to do
        otherwise
            error('filt_car_: invalid filter option (vcFilter=''%s'')', P.vcFilter);
    end % switch

    % trim padding
    if fTrim_pad
        filteredTraces = filteredTraces(n_pre + 1:end - n_post, :);
    end

    if ~getOr(P, 'fImportKilosort', 0) % global subtraction before
        [filteredTraces, vnWav2_mean] = wav_car_(filteredTraces, P);
    else
        vnWav2_mean = []; % not used when importing KiloSort data
    end
end % function
