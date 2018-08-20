%--------------------------------------------------------------------------
% 01/08/18: read and reduce
function mn = load_bin_paged_(P, viChan, nBytes_page)
    % ~35x slower than RAM indexing
    % mn = load_bin_reduce_(P, 0)
    % mn = load_bin_reduce_(P, viChans)
    % mn = load_bin_reduce_(P, viChans)
    LOAD_FACTOR = 5;
    if nargin<3, nBytes_page = []; end
    if isempty(nBytes_page)
        S = memory();
        nBytes_page = floor(S.MaxPossibleArrayBytes() / LOAD_FACTOR);
    end
    bytesPerSample = bytesPerSample_(P.dataType);
    nSamples_page = floor(nBytes_page / P.nChans / bytesPerSample);
    mn = [];

    % Determine number of samples
    if ~fileExists(P.vcFile) || isempty(viChan), return; end
    nBytes = getBytes_(P.vcFile);
    if isempty(nBytes), return; end
    headerOffset = get_(P, 'headerOffset', 0);
    nSamples = floor((nBytes-headerOffset) / bytesPerSample / P.nChans);

    % Loading loop
    fid = fopen(P.vcFile, 'r');
    try
        if headerOffset>0, fseek(fid, headerOffset, 'bof'); end
        nPages = ceil(nSamples / nSamples_page);
        if viChan(1) == 0
            fMean = 1;
            viChan = P.chanMap;
            viChan(P.viSiteZero) = [];
        else
            fMean = 0;
        end
        if nPages == 1
            mn = fread(fid, [P.nChans, nSamples], ['*', lower(P.dataType)]);
            mn = mn(viChan,:);
            if fMean, mn = cast(mean(mn), P.dataType); end
            mn = mn';
        else
            if fMean
                mn = zeros([nSamples, 1], P.dataType);
            else
                mn = zeros([nSamples, numel(viChan)], P.dataType);
            end
            for iPage = 1:nPages
                if iPage < nPages
                    nSamples_ = nSamples_page;
                else
                    nSamples_ = nSamples - nSamples_page * (iPage-1);
                end
                vi_ = (1:nSamples_) + (iPage-1) * nSamples_page;
                mn_ = fread(fid, [P.nChans, nSamples_], ['*', lower(P.dataType)]);
                mn_ = mn_(viChan,:);
                if fMean, mn_ = cast(mean(mn_), P.dataType); end
                mn(vi_,:) = mn_';
            end
        end
    catch
        disperr_();
    end
    fclose_(fid);
end %func
