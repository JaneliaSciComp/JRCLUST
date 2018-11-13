classdef Recording < handle
    %RECORDING Model of a single recording
    properties (SetAccess=private, SetObservable, Hidden, Transient)
        isError;        
        isOpen;
    end

    properties (SetAccess=private, SetObservable, Transient)
        rawData;
        filteredData;
    end

    properties (SetAccess=private, SetObservable)
        binpath;      % absolute path to binary file
        metapath;     % absolute path to meta file, if there is one

        startTime;    % beginning of recording, in samples
        endTime;      % end of recording, in samples

        dtype;        % data type contained in file
        dshape;       % shape of data (rows x columns), in samples
        headerOffset; % number of bytes at the beginning of the file to skip
    end

    properties
        spikeTimes;   % times of detected spikes in this file, in samples
        spikeSites;   % sites of detected spikes in this file, in samples
    end

    % LIFECYCLE
    methods
        function obj = Recording(filename, dtype, dshape, headerOffset)
            %RECORDING Construct an instance of this class
            % check filename exists
            obj.binpath = jrclust.utils.absPath(filename); % returns empty if not found
            obj.isError = isempty(obj.binpath);
            if obj.isError
                return;
            end

            % set object data type
            legalTypes = {'int16', 'uint16', 'int32', 'uint32', 'single', 'double'};
            if ~sum(strcmp(dtype, legalTypes)) == 1
                obj.isError = true;
                return;
            end
            obj.dtype = dtype;

            % set headerOffset
            assert(jrclust.utils.isscalarnum(headerOffset) && headerOffset >= 0, 'headerOffset must be nonnegative scalar');
            obj.headerOffset = headerOffset;

            % set object data shape
            d = dir(obj.binpath);
            sizePred = jrclust.utils.ismatrixnum(dshape) && ~isempty(dshape) ...
                && all(size(dshape) == [1 2]) && all(dshape > 0) ...
                && prod(dshape)*jrclust.utils.typeBytes(obj.dtype) == d.bytes - obj.headerOffset;
            
            if ~sizePred
                obj.isError = true;
                return;
            end
            obj.dshape = dshape;

            % load start and end times
            obj.metapath = jrclust.utils.absPath(strrep(filename, '.bin', '.meta'));
            if ~isempty(obj.metapath)
                md = jrclust.utils.metaToStruct(obj.metapath);
                obj.startTime = md.firstSample;
                obj.endTime = md.firstSample + obj.dshape(2);
            end

            obj.isOpen = false;
        end
    end

    % USER METHODS
    methods
        function open(obj)
            %OPEN Open the file for reading
            if obj.isOpen
                return;
            end

            obj.rawData = memmapfile(obj.binpath, 'Offset', obj.headerOffset, ...
                                     'Format', {obj.dtype, obj.dshape, 'Data'}, ...
                                     'Writable', false);
            obj.isOpen = true;
        end

        function close(obj)
            %CLOSE Close the file, clear its data
            if ~obj.isOpen
                return;
            end
            obj.rawData = [];
            obj.isOpen = false;
        end

        function roi = readROI(obj, rows, cols)
            %READROI get a region of interest by rows/cols
            if jrclust.utils.isscalarnum(rows)
                rows = rows*[1 1];
            end

            if jrclust.utils.isscalarnum(cols)
                cols = cols*[1 1];
            end

            rowPred = jrclust.utils.ismatrixnum(rows) && ~isempty(rows) ...
                && all(rows > 0) && issorted(rows) ...
                && rows(end) <= obj.dshape(1);
            
            colPred = jrclust.utils.ismatrixnum(cols) && ~isempty(cols) ...
                && all(cols > 0) && issorted(cols) ...
                && cols(end) <= obj.dshape(2);

            assert(rowPred, 'malformed rows');
            assert(colPred, 'malformed cols');

            if obj.isOpen % already open, don't close after read
                doClose = false;
            else
                obj.open();
                doClose = true;
            end

            roi = obj.rawData(rows, cols);

            if doClose
                obj.close();
            end
        end
    end

    % UTILITY METHODS
    methods (Access=protected, Hidden)
        function nBytes_load = file_trim_(obj, nBytes_load, P) % loadTimeLimits) % TODO: pass loadTimeLimits into this, in samples, from DetectionController
            nSamples = obj.dshape(2);

            % Apply limit to the range of samples to load
            nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
            nSamples_load = diff(nlim_load) + 1;

            nBytes_load = nSamples_load * bytesPerSample * P.nChans;
            % if nlim_load(1)>1,
            fseek_(fid, nlim_load(1), P);
            % end
        end %func

    end

    % GETTERS/SETTERS
    methods
        function data = get.rawData(obj)
            if obj.isOpen
                data = obj.rawData.Data.Data;
            else
                data = [];
            end 
        end
    end
end

