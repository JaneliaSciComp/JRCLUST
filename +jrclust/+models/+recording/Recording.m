classdef Recording < handle & dynamicprops
    %RECORDING Model of a single recording

    properties (Hidden, SetAccess=protected, SetObservable)
        hCfg;           % Config object
    end

    properties (SetAccess=private, SetObservable, Hidden, Transient)
        errMsg;
        isError;
        isOpen;
    end

    properties (SetAccess=private, SetObservable, Transient)
        rawData;      % raw data in memmapped form
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

        dataType;        % data type contained in file
        dshape;       % shape of data (rows x columns), in samples
        fSizeBytes;   % size of the file, in bytes
        headerOffset; % number of bytes at the beginning of the file to skip
    end

    properties
        spikeTimes;   % times of detected spikes in this file, in samples (relative to beginning of file)
        spikeSites;   % sites of detected spikes in this file, in samples (relative to beginning of file)
    end

    %% LIFECYCLE
    methods
        function obj = Recording(filename, hCfg) % dataType, nChans, headerOffset, 
            %RECORDING Construct an instance of this class
            % check filename exists
            obj.binpath = jrclust.utils.absPath(filename); % returns empty if not found
            obj.isError = isempty(obj.binpath);
            if obj.isError
                obj.errMsg = sprintf('file not found: %s', filename);
                return;
            end

            % set object data type
            obj.dataType = hCfg.dataType;

            % set headerOffset
            obj.headerOffset = hCfg.headerOffset;

            % set object data shape
            d = dir(obj.binpath);
            obj.fSizeBytes = d.bytes;

            nSamples = (d.bytes - obj.headerOffset) / hCfg.nChans / jrclust.utils.typeBytes(obj.dataType);
            if ceil(nSamples) ~= nSamples % must be an integer or we don't have a complete recording
                obj.errMsg = 'Number of samples computed is not an integer. Check your sample rate or nChans?';
                obj.isError = 1;
                return;
            end
            obj.dshape = [hCfg.nChans nSamples];

            % load start and end times
            obj.metapath = jrclust.utils.findMeta(filename);
            if ~isempty(obj.metapath)
                md = jrclust.utils.metaToStruct(obj.metapath);
                if isfield(md, 'firstSample')
                    obj.startTime = md.firstSample;
                    obj.endTime = md.firstSample + obj.dshape(2);
                end
            end

            try
                obj.hCfg = hCfg;
            catch ME
                obj.errMsg = ME.message;
                obj.isError = 1;
                return;
            end

            obj.isOpen = 0;
        end
    end

    %% USER METHODS
    methods
        function open(obj)
            %OPEN Open the file for reading
            if obj.isOpen
                return;
            end

            obj.rawData = memmapfile(obj.binpath, 'Offset', obj.headerOffset, ...
                                     'Format', {obj.dataType, obj.dshape, 'Data'}, ...
                                     'Writable', false);
            obj.isOpen = 1;
        end

        function close(obj)
            %CLOSE Close the file, clear its data
            if ~obj.isOpen
                return;
            end
            obj.rawData = [];
            obj.isOpen = 0;
        end

        function roi = readROI(obj, rows, cols)
            %READROI Get a region of interest by rows/cols
            rowPred = jrclust.utils.ismatrixnum(rows) && ~isempty(rows) ...
                && all(rows > 0) && max(rows) <= obj.dshape(1);
            
            colPred = jrclust.utils.ismatrixnum(cols) && ~isempty(cols) ...
                && all(cols > 0) && issorted(cols) ...
                && cols(end) <= obj.dshape(2);

            assert(rowPred, 'malformed rows');
            assert(colPred, 'malformed cols');

            if obj.isOpen % already open, don't close after read
                doClose = 0;
            else
                obj.open();
                doClose = 1;
            end

            roi = obj.rawData(rows, cols);

            if doClose
                obj.close();
            end
        end

        function setDetections(obj, spikeData)
            %SETDETECTIONS Assign spikeTimes and spikeSites, optionally other fields
            % reject spikeData without spikeTimes or spikeSites
            if ~isstruct(spikeData) || ~isfield(spikeData, 'spikeTimes') || ~isfield(spikeData, 'spikeSites')
                error('spikeData must be a struct containing, at minimum, spikeTimes and spikeSites');
            end

            % reject spikeData if counts mismatch
            if numel(spikeData.spikeTimes) ~= numel(spikeData.spikeSites)
                error('counts of spikeTimes and spikeSites mismatch');
            end

            % reject spikeData if spikeTimes or spikeSites go out of range
            if any(spikeData.spikeTimes < 1) || any(spikeData.spikeTimes > obj.dshape(2))
                error('spikeTimes out of range');
            end
            if any(spikeData.spikeSites < 1) || any(spikeData.spikeSites > obj.dshape(1))
                error('spikeSites out of range');
            end

            obj.spikeTimes = spikeData.spikeTimes;
            spikeData = rmfield(spikeData, 'spikeTimes');

            obj.spikeSites = spikeData.spikeSites;
            spikeData = rmfield(spikeData, 'spikeSites');

            if ~isempty(spikeData)
                fieldNames = fieldnames(spikeData);
                for i = 1:numel(fieldNames)
                    fn = fieldNames{i};
                    if ismember(fn, {'spikesRaw', 'spikesFilt', 'spikeFeatures'})
                        continue; % don't save these twice
                    end

                    if ~isprop(obj, fn) % don't overwrite any fields already here
                        obj.addprop(fn);
                        obj.(fn) = spikeData.(fn);
                    end
                end
            end
        end
    end

    %% GETTERS/SETTERS
    methods
        % hCfg
        function set.hCfg(obj, hc)
            failMsg = 'hCfg must be an instance of jrclust.Config';
            assert(isempty(hc) || isa(hc, 'jrclust.Config'), failMsg);
            obj.hCfg = hc;
        end

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

