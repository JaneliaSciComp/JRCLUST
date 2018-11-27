classdef Recording < handle
    %RECORDING Model of a single recording
    properties (SetAccess=private, SetObservable, Hidden, Transient)
        errMsg;
        isError;        
        isOpen;
    end

    properties (SetAccess=private, SetObservable, Transient)
        rawData;
        filteredData;
    end

    properties (SetObservable, Dependent, Transient)
        nChans;       % number of channels (rows) in the recording
        nSamples;     % number of samples per channel (columns) in the recording
    end

    properties (SetAccess=private, SetObservable)
        binpath;      % absolute path to binary file
        metapath;     % absolute path to meta file, if there is one

        startTime;    % beginning of recording, in samples
        endTime;      % end of recording, in samples

        dtype;        % data type contained in file
        dshape;       % shape of data (rows x columns), in samples
        fSizeBytes;   % size of the file, in bytes
        headerOffset; % number of bytes at the beginning of the file to skip
    end

    properties
        spikeTimes;   % times of detected spikes in this file, in samples
        spikeSites;   % sites of detected spikes in this file, in samples
    end

    % LIFECYCLE
    methods
        function obj = Recording(filename, dtype, nChans, headerOffset)
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
                obj.errMsg = sprintf('dtype must be one of %s', strjoin(legalTypes, ', '));
                obj.isError = true;
                return;
            end
            obj.dtype = dtype;

            % set headerOffset
            if ~(jrclust.utils.isscalarnum(headerOffset) && headerOffset >= 0)
                obj.errMsg = 'headerOffset must be nonnegative scalar';
                obj.isError = true;
                return;
            end
            obj.headerOffset = headerOffset;

            % set object data shape
            if ~(jrclust.utils.isscalarnum(nChans) && nChans > 0)
                obj.errMsg = 'nChans must be a positive scalar';
                obj.isError = true;
                return;
            end
            d = dir(obj.binpath);
            obj.fSizeBytes = d.bytes;

            nSamples = (d.bytes - obj.headerOffset) / nChans / jrclust.utils.typeBytes(obj.dtype);
            if ceil(nSamples) ~= nSamples % must be an integer or we don't have a complete recording
                obj.errMsg = 'incomplete recording';
                obj.isError = true;
                return;
            end
            obj.dshape = [nChans nSamples];

            % load start and end times
            obj.metapath = jrclust.utils.findMeta(filename);
            if ~isempty(obj.metapath)
                md = jrclust.utils.metaToStruct(obj.metapath);
                if isfield(md, 'firstSample')
                    obj.startTime = md.firstSample;
                    obj.endTime = md.firstSample + obj.dshape(2);
                end
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
                && all(rows > 0) && max(rows) <= obj.dshape(1);
            
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
    methods (Hidden)
        function nBytesLoad = subsetBytes(obj, loadTimeLimits)
            %SUBSETBYTES Get number of bytes to load from file, given time limits
            nBytesLoad = obj.fSizeBytes - obj.headerOffset;

            if isempty(loadTimeLimits)
                return;
            end

            loadLimits = min(max(loadTimeLimits, 1), obj.nSamples);
            nSamplesLoad = diff(loadLimits) + 1;
            nBytesLoad = nSamplesLoad * jrclust.utils.typeBytes(obj.dtype) * obj.nChans;
        end
    end

    % GETTERS/SETTERS
    methods
        % rawData
        function rd = get.rawData(obj)
            if obj.isOpen
                rd = obj.rawData.Data.Data;
            else
                rd = [];
            end
        end

        % nChans
        function nc = get.nChans(obj)
            nc = obj.dshape(1);
        end

        % nSamples
        function nc = get.nSamples(obj)
            nc = obj.dshape(2);
        end
    end
end

