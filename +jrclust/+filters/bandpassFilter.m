function samplesOut = bandpassFilter(samplesIn, hCfg)
    %BANDPASSFILTER Summary of this function goes here
    try
        hCfg.updateLog('bpFilt', sprintf('Applying bandpass filter to %d samples', size(samplesIn, 1)), 1, 0);
        samplesOut = filtfiltChain(single(samplesIn), hCfg);
    catch ME
        hCfg.updateLog('bpFilt', sprintf('GPU filtering failed: %s (retrying on CPU)', ME.message), 1, 0);
        hCfg.useGPUFilt = 0;
        samplesOut = filtfiltChain(single(samplesIn), hCfg);
    end

    hCfg.updateLog('bpFilt', 'Finished filtering', 0, 1);

    samplesOut = int16(samplesOut);
end

%% LOCAL FUNCTIONS
function samplesIn = filtfiltChain(samplesIn, filtOpts)
    %FILTFILTCHAIN Construct a filter chain
    [cvrA, cvrB] = deal({});

    if ~isempty(filtOpts.freqLimBP)
       [cvrB{end+1}, cvrA{end+1}] = makeFilt_(filtOpts.freqLimBP, 'bandpass', filtOpts);
    end

    if ~isempty(filtOpts.freqLimStop)
       [cvrB{end+1}, cvrA{end+1}] = makeFilt_(filtOpts.freqLimStop, 'stop', filtOpts);
    end  

    if ~isempty(filtOpts.freqLimNotch)
        if ~iscell(filtOpts.freqLimNotch)
            [cvrB{end+1}, cvrA{end+1}] = makeFilt_(filtOpts.freqLimNotch, 'notch', filtOpts);
        else
            for iCell=1:numel(filtOpts.freqLimNotch)
                [cvrB{end+1}, cvrA{end+1}] = makeFilt_(filtOpts.freqLimNotch{iCell}, 'notch', filtOpts);
            end
        end
    end    

    %----------------
    % Run the filter chain
    fInt16 = isa(samplesIn, 'int16');
    if filtOpts.useGPUFilt
        samplesIn = gpuArray(samplesIn);
    end

    samplesIn = filt_pad_('add', samplesIn, filtOpts.nSamplesPad); % slow
    if fInt16
        samplesIn = single(samplesIn); % double for long data?
    end

    % first pass
    for iFilt=1:numel(cvrA)
        samplesIn = flipud(filter(cvrB{iFilt}, cvrA{iFilt}, ...
                flipud(filter(cvrB{iFilt}, cvrA{iFilt}, samplesIn))));
    end

    if filtOpts.gainBoost ~= 1
        samplesIn = samplesIn * filtOpts.gainBoost;
    end

    if fInt16
        samplesIn = int16(samplesIn);
    end

    samplesIn = filt_pad_('remove', samplesIn, filtOpts.nSamplesPad); % slow    
end

function [vrFiltB, vrFiltA] = makeFilt_(freqLimBP, vcType, filtOpts)
    if nargin < 2
        vcType = 'bandpass';
    end

    freqLimBP = freqLimBP / filtOpts.sampleRate * 2;

    if ~strcmpi(vcType, 'notch')
        if filtOpts.useElliptic  % copied from wave_clus
            if isinf(freqLimBP(1)) || freqLimBP(1) <= 0
                [vrFiltB, vrFiltA]=ellip(filtOpts.filtOrder,0.1,40, freqLimBP(2), 'low');
            elseif isinf(freqLimBP(2))
                [vrFiltB, vrFiltA]=ellip(filtOpts.filtOrder,0.1,40, freqLimBP(1), 'high');
            else
                [vrFiltB, vrFiltA]=ellip(filtOpts.filtOrder,0.1,40, freqLimBP, vcType);
            end    
        else
            if isinf(freqLimBP(1)) || freqLimBP(1) <= 0
                [vrFiltB, vrFiltA] = butter(filtOpts.filtOrder, freqLimBP(2),'low');        
            elseif isinf(freqLimBP(2))
                [vrFiltB, vrFiltA] = butter(filtOpts.filtOrder, freqLimBP(1),'high');
            else
                [vrFiltB, vrFiltA] = butter(filtOpts.filtOrder, freqLimBP, vcType);    
            end
        end
    else
        [vrFiltB, vrFiltA] = iirNotch(mean(freqLimBP), diff(freqLimBP));
    end

    if filtOpts.useGPUFilt  
        vrFiltB = gpuArray(vrFiltB);
        vrFiltA = gpuArray(vrFiltA);
    end
end

function mrWav = filt_pad_(vcMode, mrWav, nPad)
    % add padding by reflection. removes artifact at the end

    if isempty(nPad), return; end
    if nPad==0, return; end
    nPad = min(nPad, size(mrWav,1));

    switch lower(vcMode)
        case 'add'
            mrWav = [flipud(mrWav(1:nPad,:)); mrWav; flipud(mrWav(end-nPad+1:end,:))];
        case 'remove'
            mrWav = mrWav(nPad+1:end-nPad,:);
    end %switch
end

function [num,den] = iirNotch(Wo, BW)
    % Define default values.
    Ab = abs(10*log10(.5)); % 3-dB width
    % Design a 2nd-order notch digital filter.

    % Inputs are normalized by pi.
    BW = BW*pi;
    Wo = Wo*pi;

    Gb   = 10^(-Ab/20);
    beta = (sqrt(1-Gb.^2)/Gb)*tan(BW/2);
    gain = 1/(1+beta);

    num  = gain*[1 -2*cos(Wo) 1];
    den  = [1 -2*gain*cos(Wo) (2*gain-1)];
end