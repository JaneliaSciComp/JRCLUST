%--------------------------------------------------------------------------
function [vrPow, vrFreq] = plotMedPower_(mrData, varargin)

    if numel(varargin) ==1
        if ~isstruct(varargin{1}), P.sampleRateHz = varargin{1};
            else P = varargin{1};
        end
    else
        P = funcInStr_(varargin{:});
    end

    P = funcDefStr_(P, 'viChanExcl', [1 18 33 50 65 82 97 114], 'sampleRateHz', 25000, ...
    'nSmooth', 3, 'LineStyle', 'k-', 'fPlot', 1, 'fKHz', 0, 'vcMode', 'max');

    if iscell(mrData) %batch mode
        csFname = mrData;
        [vrPow, vrFreq] =deal(cell(numel(csFname)));
        for i=1:numel(csFname)
            hold on;
            [vrPow{i}, vrFreq{i}] = plotMedPower_(csFname{i});
        end
        legend(csFname);
        return;
    else
        vcFname='';
    end
    warning off;

    mrData = fft(mrData);
    mrData = real(mrData .* conj(mrData)) / size(mrData,1) / (P.sampleRateHz/2);
    % mrPowFilt = filter([1 1 1], 3, mrPow);
    % mrPow = fftshift(mrPow);
    imid0 = ceil(size(mrData,1)/2);
    vrFreq = (0:size(mrData,1)-1) * (P.sampleRateHz/size(mrData,1));
    vrFreq = vrFreq(2:imid0);
    viChan = setdiff(1:size(mrData,2), P.viChanExcl);
    if size(mrData,2)>1
        switch P.vcMode
            case 'mean'
            vrPow = mean(mrData(2:imid0,viChan), 2);
            case 'max'
            vrPow = max(mrData(2:imid0,viChan), [], 2);
        end
    else
        vrPow = mrData(2:imid0,1);
    end
    % vrPow = std(mrData(:,viChan), 1, 2);
    if P.nSmooth>1, vrPow = filterq_(ones([P.nSmooth,1]),P.nSmooth,vrPow); end

    if P.fPlot
        if P.fKHz, vrFreq = vrFreq/1000; end
        plot(vrFreq, pow2db_(vrPow), P.LineStyle);
        %     set(gca, 'YScale', 'log');
        %     set(gca, 'XScale', 'linear');
        xlabel('Frequency (Hz)'); ylabel('Mean power across sites (dB uV^2/Hz)');
        % xlim_([0 P.sampleRateHz/2]);
        grid on;
        try
            xlim_(vrFreq([1, end]));
            set(gcf,'color','w');
            title(vcFname, 'Interpreter', 'none');
        catch

        end
    end
end %func
