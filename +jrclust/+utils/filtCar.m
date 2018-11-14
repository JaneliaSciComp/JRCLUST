%--------------------------------------------------------------------------
% function [samplesIn, vnWav2_mean] = filtCar(samplesIn, hCfg, windowPre, windowPost, fTrim_pad)
function [samplesIn, vnWav2_mean] = filtCar(samplesIn, windowPre, windowPost, trimPad, hCfg)
    %FILTCAR Apply user-specified filter and common-average referencing

    if nargin < 3
        windowPre = [];
    end
    if nargin < 4
        windowPost = [];
    end
    if nargin < 5
        trimPad = true;
    end

    nPadPre = size(windowPre, 1);
    nPadPost = size(windowPost, 1);

    samplesIn = [windowPre; samplesIn; windowPost];

    % apply filter
    if strcmp(hCfg.filterType, 'user') && ~isempty(hCfg.userFiltKernel)
        for i = 1:size(samplesIn, 2)
            samplesIn(:, i) = conv(samplesIn(:,i), hCfg.userFiltKernel, 'same');
        end
    elseif strcmp(hCfg.filterType, 'fir1')
        n5ms = round(hCfg.sRateHz / 1000 * 5);

        firFilt = single(fir1(n5ms, hCfg.freqLim/hCfg.sRateHz*2));
        for i = 1:size(samplesIn, 2)
            samplesIn(:,i) = conv(samplesIn(:, i), firFilt, 'same');
        end
    elseif strcmp(hCfg.filterType, 'ndiff')
        samplesIn = ndiff(samplesIn, hCfg.nDiff_filt);
    end

%         case 'ndiff'
%             samplesIn = ndiff(samplesIn, hCfg.nDiff_filt);
% 
%         case 'fftdiff'
%             samplesIn = fftdiff_(samplesIn, hCfg);
% 
%         case {'sgdiff', 'sgfilt'}
%             samplesIn = sgfilt_(samplesIn, hCfg.nDiff_filt);
% 
%         case 'bandpass'
%             try
%                 samplesIn = filtfilt_chain(single(samplesIn), hCfg);
%             catch
%                 fprintf('GPU filtering failed. Trying CPU filtering.\n');
%                 samplesIn = filtfilt_chain(single(samplesIn), setfield(hCfg, 'fGpu', 0));
%             end
%             samplesIn = int16(samplesIn);
% 
%         case {'none', 'skip'} % no filter is applied
%             
%         case 'ndist'
%             samplesIn = ndist_filt_(samplesIn, get_set_(hCfg, 'ndist_filt', 5));
% 
%         otherwise
%             error('invalid filter option (filterType=''%s'')', filterType);
%     end % switch

    % trim padding
    if (nPadPre > 0 || nPadPost > 0) && trimPad
        samplesIn = samplesIn(nPadPre+1:end-nPadPost,:);
    end

    % global subtraction before
    [samplesIn, vnWav2_mean] = wav_car_(samplesIn, hCfg);
end

%% LOCAL FUNCTIONS
function mnWav = ndiff(mnWav, nDiff_filt)
    % works for a vector, matrix and tensor. Like sgdiff but averaging instead
    % fInvert_filter = 0;

    % n1 = size(mnWav,1);
    % if n1==1, n1 = size(mnWav,2);  end
    % if nDiff_filt==0, mnWav1 = mnWav; return; end
    % if fInvert_filter
    %     [miB, miA] = sgfilt_init_(n1, nDiff_filt, isGpu_(mnWav));
    % else
    %     [miA, miB] = sgfilt_init_(n1, nDiff_filt, isGpu_(mnWav));
    % end
    if isempty(intersect(nDiff_filt, [1 2 3]))
        return;
    end

    if nDiff_filt == 1
        mnWav(1:end-1,:) = diff(mnWav);
        mnWav(end,:) = 0;
    elseif nDiff_filt == 2
        mnWav(2:end-2,:) = 2*diff(mnWav(2:end-1,:)) + (mnWav(4:end,:) - mnWav(1:end-3,:));
        mnWav([1,end-1,end],:) = 0;
    else % nDiff_filt == 3
        mnWav(3:end-3,:) = 3*diff(mnWav(3:end-2,:)) + 2*(mnWav(5:end-1,:) - mnWav(2:end-4,:)) + (mnWav(6:end,:) - mnWav(1:end-5,:));
        mnWav([1,2,end-2,end-1,end],:) = 0;
    end
end
