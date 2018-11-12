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

        startTime;    % beginning of recording, in seconds
        endTime;      % end of recording, in seconds

        dtype;        % data type contained in file
        dshape;       % shape of data (rows x columns), in samples
        headerOffset; % number of bytes at the beginning of the file to skip
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

            % TODO: if meta file exists, load up metadata about it
            obj.metapath = jrclust.utils.absPath(strrep(filename, '.bin', '.meta'));

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
            if ~obj.isOpen
                return;
            end
            obj.rawData = [];
            obj.isOpen = false;
        end

        function roi = readROI(obj, rowBounds, colBounds)
            if jrclust.utils.isscalarnum(rowBounds)
                rowBounds = rowBounds*[1 1];
            end

            if jrclust.utils.isscalarnum(colBounds)
                colBounds = colBounds*[1 1];
            end

            rowPred = jrclust.utils.ismatrixnum(rowBounds) && ~isempty(rowBounds) ...
                && all(size(rowBounds) == [1 2]) && all(rowBounds > 0) ...
                && rowBounds(2) <= obj.dshape(1);
            
            colPred = jrclust.utils.ismatrixnum(colBounds) && ~isempty(colBounds) ...
                && all(size(colBounds) == [1 2]) && all(colBounds > 0) ...
                && colBounds(2) <= obj.dshape(2);

            assert(rowPred, 'malformed rowBounds');
            assert(colPred, 'malformed colBounds');

            if obj.isOpen % already open, don't close after read
                doClose = false;
            else
                obj.open();
                doClose = true;
            end

            roi = obj.rawData(rowBounds(1):rowBounds(2), colBounds(1):colBounds(2));
            
            if doClose
                obj.close();
            end
        end
    end

    % GETTERS/SETTERS
    methods
        function data = get.rawData(obj)
            data = obj.rawData.Data.Data;
        end
    end
end

