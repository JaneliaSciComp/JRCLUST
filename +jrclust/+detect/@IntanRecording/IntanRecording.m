classdef IntanRecording < jrclust.interfaces.RawRecording
    %INTANRECORDING Model of a single recording in the traditional Intan file format
    %% INTAN-SPECIFIC PROPERTIES
    properties (Hidden, SetAccess=private)
        rawFid;
    end

    properties (Hidden, SetObservable)
        fileData;
    end

    properties (Dependent, SetAccess=private, Transient)
        ampData;
        ampChannels;
        auxInputChannels;
        auxInputData;

        boardADCChannels;
        boardADCData;
        boardDigInChannels;
        boardDigInData;
        boardDigInRaw;
        boardDigOutChannels;
        boardDigOutData;
        boardDigOutRaw;
        boardMode;

        dspCutoffFreqActual;
        dspCutoffFreqRequested;
        dspEnabled;

        impedanceTestFreqActual;
        impedanceTestFreqRequested;

        lowerBandwidthActual;
        lowerBandwidthRequested;

        nAmpChannels;
        nAmpSamples;
        nAuxInputChannels;
        nBoardADCChannels;
        nBoardDigInChannels;
        nBoardDigOutChannels;
        nBytesBlock;
        nSamplesBlock;
        nSupplyVoltageChannels;
        nTempSensors;
        notchFilterFreq;
        notchFilterMode;
        notes;

        refChannelName;

        sampleRate;
        supplyVoltageChannels;
        supplyVoltageData;

        tAmpData;
        tAuxInput;
        tBoardADC;
        tDig;
        tSupplyVoltage;
        tTempSensor;
        tempSensorData;

        upperBandwidthActual;
        upperBandwidthRequested;

        versionMajor;
        versionMinor;
    end

    %% LIFECYCLE    
    methods
        function obj = IntanRecording(filename, hCfg)
            %INTANRECORDING Construct an instance of this class
            % set object data type
            obj = obj@jrclust.interfaces.RawRecording(filename, hCfg);
            obj.dataType = hCfg.dataType;

            try
                obj.hCfg = hCfg;
            catch ME
                obj.errMsg = ME.message;
                obj.isError = 1;
                return;
            end

            % set headerOffset, dshape
            obj.fileData = obj.readHeader(obj.rawPath);
            obj.headerOffset = obj.fileData.headerOffset;
            obj.dshape = [obj.fileData.nChans, obj.fileData.nAmpSamples];

            obj.rawFid = -1;
            obj.rawIsOpen = 0;
        end
    end

    %% USER METHODS
    methods (Static)
        fileData = readHeader(rawPath);
    end

    %% GETTERS/SETTERS
    methods
        function val = get.ampData(obj)
            if isfield(obj.fileData, 'ampData')
                val = obj.fileData.ampData;
            else
                val = [];
            end
        end
        function set.ampData(obj, val)
            obj.fileData.ampData = val;
        end

        function val = get.ampChannels(obj)
            if isfield(obj.fileData, 'ampChannels')
                val = obj.fileData.ampChannels;
            else
                val = [];
            end
        end
        function set.ampChannels(obj, val)
            obj.fileData.ampChannels = val;
        end

        function val = get.auxInputChannels(obj)
            if isfield(obj.fileData, 'auxInputChannels')
                val = obj.fileData.auxInputChannels;
            else
                val = [];
            end
        end
        function set.auxInputChannels(obj, val)
            obj.fileData.auxInputChannels = val;
        end

        function val = get.auxInputData(obj)
            if isfield(obj.fileData, 'auxInputData')
                val = obj.fileData.auxInputData;
            else
                val = [];
            end
        end
        function set.auxInputData(obj, val)
            obj.fileData.auxInputData = val;
        end

        function val = get.boardADCChannels(obj)
            if isfield(obj.fileData, 'boardADCChannels')
                val = obj.fileData.boardADCChannels;
            else
                val = [];
            end
        end
        function set.boardADCChannels(obj, val)
            obj.fileData.boardADCChannels = val;
        end

        function val = get.boardADCData(obj)
            if isfield(obj.fileData, 'boardADCData')
                val = obj.fileData.boardADCData;
            else
                val = [];
            end
        end
        function set.boardADCData(obj, val)
            obj.fileData.boardADCData = val;
        end

        function val = get.boardDigInChannels(obj)
            if isfield(obj.fileData, 'boardDigInChannels')
                val = obj.fileData.boardDigInChannels;
            else
                val = [];
            end
        end
        function set.boardDigInChannels(obj, val)
            obj.fileData.boardDigInChannels = val;
        end

        function val = get.boardDigInData(obj)
            if isfield(obj.fileData, 'boardDigInData')
                val = obj.fileData.boardDigInData;
            else
                val = [];
            end
        end
        function set.boardDigInData(obj, val)
            obj.fileData.boardDigInData = val;
        end

        function val = get.boardDigInRaw(obj)
            if isfield(obj.fileData, 'boardDigInRaw')
                val = obj.fileData.boardDigInRaw;
            else
                val = [];
            end
        end
        function set.boardDigInRaw(obj, val)
            obj.fileData.boardDigInRaw = val;
        end

        function val = get.boardDigOutChannels(obj)
            if isfield(obj.fileData, 'boardDigOutChannels')
                val = obj.fileData.boardDigOutChannels;
            else
                val = [];
            end
        end
        function set.boardDigOutChannels(obj, val)
            obj.fileData.boardDigOutChannels = val;
        end

        function val = get.boardDigOutData(obj)
            if isfield(obj.fileData, 'boardDigOutData')
                val = obj.fileData.boardDigOutData;
            else
                val = [];
            end
        end
        function set.boardDigOutData(obj, val)
            obj.fileData.boardDigOutData = val;
        end

        function val = get.boardDigOutRaw(obj)
            if isfield(obj.fileData, 'boardDigOutRaw')
                val = obj.fileData.boardDigOutRaw;
            else
                val = [];
            end
        end
        function set.boardDigOutRaw(obj, val)
            obj.fileData.boardDigOutRaw = val;
        end

        function val = get.boardMode(obj)
            if isfield(obj.fileData, 'boardMode')
                val = obj.fileData.boardMode;
            else
                val = [];
            end
        end
        function set.boardMode(obj, val)
            obj.fileData.boardMode = val;
        end

        function val = get.dspCutoffFreqActual(obj)
            if isfield(obj.fileData, 'dspCutoffFreqActual')
                val = obj.fileData.dspCutoffFreqActual;
            else
                val = [];
            end
        end
        function set.dspCutoffFreqActual(obj, val)
            obj.fileData.dspCutoffFreqActual = val;
        end

        function val = get.dspCutoffFreqRequested(obj)
            if isfield(obj.fileData, 'dspCutoffFreqRequested')
                val = obj.fileData.dspCutoffFreqRequested;
            else
                val = [];
            end
        end
        function set.dspCutoffFreqRequested(obj, val)
            obj.fileData.dspCutoffFreqRequested = val;
        end

        function val = get.dspEnabled(obj)
            if isfield(obj.fileData, 'dspEnabled')
                val = obj.fileData.dspEnabled;
            else
                val = [];
            end
        end
        function set.dspEnabled(obj, val)
            obj.fileData.dspEnabled = val;
        end

        function val = get.impedanceTestFreqActual(obj)
            if isfield(obj.fileData, 'impedanceTestFreqActual')
                val = obj.fileData.impedanceTestFreqActual;
            else
                val = [];
            end
        end
        function set.impedanceTestFreqActual(obj, val)
            obj.fileData.impedanceTestFreqActual = val;
        end

        function val = get.impedanceTestFreqRequested(obj)
            if isfield(obj.fileData, 'impedanceTestFreqRequested')
                val = obj.fileData.impedanceTestFreqRequested;
            else
                val = [];
            end
        end
        function set.impedanceTestFreqRequested(obj, val)
            obj.fileData.impedanceTestFreqRequested = val;
        end

        function val = get.lowerBandwidthActual(obj)
            if isfield(obj.fileData, 'lowerBandwidthActual')
                val = obj.fileData.lowerBandwidthActual;
            else
                val = [];
            end
        end
        function set.lowerBandwidthActual(obj, val)
            obj.fileData.lowerBandwidthActual = val;
        end

        function val = get.lowerBandwidthRequested(obj)
            if isfield(obj.fileData, 'lowerBandwidthRequested')
                val = obj.fileData.lowerBandwidthRequested;
            else
                val = [];
            end
        end
        function set.lowerBandwidthRequested(obj, val)
            obj.fileData.lowerBandwidthRequested = val;
        end

        function val = get.nAmpChannels(obj)
            if isfield(obj.fileData, 'nAmpChannels')
                val = obj.fileData.nAmpChannels;
            else
                val = [];
            end
        end
        function set.nAmpChannels(obj, val)
            obj.fileData.nAmpChannels = val;
        end

        function val = get.nAmpSamples(obj)
            if isfield(obj.fileData, 'nAmpSamples')
                val = obj.fileData.nAmpSamples;
            else
                val = [];
            end
        end
        function set.nAmpSamples(obj, val)
            obj.fileData.nAmpSamples = val;
        end

        function val = get.nAuxInputChannels(obj)
            if isfield(obj.fileData, 'nAuxInputChannels')
                val = obj.fileData.nAuxInputChannels;
            else
                val = [];
            end
        end
        function set.nAuxInputChannels(obj, val)
            obj.fileData.nAuxInputChannels = val;
        end

        function val = get.nBoardADCChannels(obj)
            if isfield(obj.fileData, 'nBoardADCChannels')
                val = obj.fileData.nBoardADCChannels;
            else
                val = [];
            end
        end
        function set.nBoardADCChannels(obj, val)
            obj.fileData.nBoardADCChannels = val;
        end

        function val = get.nBoardDigInChannels(obj)
            if isfield(obj.fileData, 'nBoardDigInChannels')
                val = obj.fileData.nBoardDigInChannels;
            else
                val = [];
            end
        end
        function set.nBoardDigInChannels(obj, val)
            obj.fileData.nBoardDigInChannels = val;
        end

        function val = get.nBoardDigOutChannels(obj)
            if isfield(obj.fileData, 'nBoardDigOutChannels')
                val = obj.fileData.nBoardDigOutChannels;
            else
                val = [];
            end
        end
        function set.nBoardDigOutChannels(obj, val)
            obj.fileData.nBoardDigOutChannels = val;
        end

        function val = get.nBytesBlock(obj)
            if isfield(obj.fileData, 'nBytesBlock')
                val = obj.fileData.nBytesBlock;
            else
                val = [];
            end
        end
        function set.nBytesBlock(obj, val)
            obj.fileData.nBytesBlock = val;
        end

        function val = get.nSamplesBlock(obj)
            if isfield(obj.fileData, 'nSamplesBlock')
                val = obj.fileData.nSamplesBlock;
            else
                val = [];
            end
        end
        function set.nSamplesBlock(obj, val)
            obj.fileData.nSamplesBlock = val;
        end

        function val = get.nSupplyVoltageChannels(obj)
            if isfield(obj.fileData, 'nSupplyVoltageChannels')
                val = obj.fileData.nSupplyVoltageChannels;
            else
                val = [];
            end
        end
        function set.nSupplyVoltageChannels(obj, val)
            obj.fileData.nSupplyVoltageChannels = val;
        end

        function val = get.nTempSensors(obj)
            if isfield(obj.fileData, 'nTempSensors')
                val = obj.fileData.nTempSensors;
            else
                val = [];
            end
        end
        function set.nTempSensors(obj, val)
            obj.fileData.nTempSensors = val;
        end

        function val = get.notchFilterFreq(obj)
            if isfield(obj.fileData, 'notchFilterFreq')
                val = obj.fileData.notchFilterFreq;
            else
                val = [];
            end
        end
        function set.notchFilterFreq(obj, val)
            obj.fileData.notchFilterFreq = val;
        end

        function val = get.notchFilterMode(obj)
            if isfield(obj.fileData, 'notchFilterMode')
                val = obj.fileData.notchFilterMode;
            else
                val = [];
            end
        end
        function set.notchFilterMode(obj, val)
            obj.fileData.notchFilterMode = val;
        end

        function val = get.notes(obj)
            if isfield(obj.fileData, 'notes')
                val = obj.fileData.notes;
            else
                val = [];
            end
        end
        function set.notes(obj, val)
            obj.fileData.notes = val;
        end

        function val = get.refChannelName(obj)
            if isfield(obj.fileData, 'refChannelName')
                val = obj.fileData.refChannelName;
            else
                val = [];
            end
        end
        function set.refChannelName(obj, val)
            obj.fileData.refChannelName = val;
        end

        function val = get.sampleRate(obj)
            if isfield(obj.fileData, 'sampleRate')
                val = obj.fileData.sampleRate;
            else
                val = [];
            end
        end
        function set.sampleRate(obj, val)
            obj.fileData.sampleRate = val;
        end

        function val = get.supplyVoltageChannels(obj)
            if isfield(obj.fileData, 'supplyVoltageChannels')
                val = obj.fileData.supplyVoltageChannels;
            else
                val = [];
            end
        end
        function set.supplyVoltageChannels(obj, val)
            obj.fileData.supplyVoltageChannels = val;
        end

        function val = get.supplyVoltageData(obj)
            if isfield(obj.fileData, 'supplyVoltageData')
                val = obj.fileData.supplyVoltageData;
            else
                val = [];
            end
        end
        function set.supplyVoltageData(obj, val)
            obj.fileData.supplyVoltageData = val;
        end

        function val = get.tAmpData(obj)
            if isfield(obj.fileData, 'tAmpData')
                val = obj.fileData.tAmpData;
            else
                val = [];
            end
        end
        function set.tAmpData(obj, val)
            obj.fileData.tAmpData = val;
        end

        function val = get.tAuxInput(obj)
            if isfield(obj.fileData, 'tAuxInput')
                val = obj.fileData.tAuxInput;
            else
                val = [];
            end
        end
        function set.tAuxInput(obj, val)
            obj.fileData.tAuxInput = val;
        end

        function val = get.tBoardADC(obj)
            if isfield(obj.fileData, 'tBoardADC')
                val = obj.fileData.tBoardADC;
            else
                val = [];
            end
        end
        function set.tBoardADC(obj, val)
            obj.fileData.tBoardADC = val;
        end

        function val = get.tDig(obj)
            if isfield(obj.fileData, 'tDig')
                val = obj.fileData.tDig;
            else
                val = [];
            end
        end
        function set.tDig(obj, val)
            obj.fileData.tDig = val;
        end

        function val = get.tSupplyVoltage(obj)
            if isfield(obj.fileData, 'tSupplyVoltage')
                val = obj.fileData.tSupplyVoltage;
            else
                val = [];
            end
        end
        function set.tSupplyVoltage(obj, val)
            obj.fileData.tSupplyVoltage = val;
        end

        function val = get.tTempSensor(obj)
            if isfield(obj.fileData, 'tTempSensor')
                val = obj.fileData.tTempSensor;
            else
                val = [];
            end
        end
        function set.tTempSensor(obj, val)
            obj.fileData.tTempSensor = val;
        end

        function val = get.tempSensorData(obj)
            if isfield(obj.fileData, 'tempSensorData')
                val = obj.fileData.tempSensorData;
            else
                val = [];
            end
        end
        function set.tempSensorData(obj, val)
            obj.fileData.tempSensorData = val;
        end

        function val = get.upperBandwidthActual(obj)
            if isfield(obj.fileData, 'upperBandwidthActual')
                val = obj.fileData.upperBandwidthActual;
            else
                val = [];
            end
        end
        function set.upperBandwidthActual(obj, val)
            obj.fileData.upperBandwidthActual = val;
        end

        function val = get.upperBandwidthRequested(obj)
            if isfield(obj.fileData, 'upperBandwidthRequested')
                val = obj.fileData.upperBandwidthRequested;
            else
                val = [];
            end
        end
        function set.upperBandwidthRequested(obj, val)
            obj.fileData.upperBandwidthRequested = val;
        end

        function val = get.versionMajor(obj)
            if isfield(obj.fileData, 'versionMajor')
                val = obj.fileData.versionMajor;
            else
                val = [];
            end
        end
        function set.versionMajor(obj, val)
            obj.fileData.versionMajor = val;
        end

        function val = get.versionMinor(obj)
            if isfield(obj.fileData, 'versionMinor')
                val = obj.fileData.versionMinor;
            else
                val = [];
            end
        end
        function set.versionMinor(obj, val)
            obj.fileData.versionMinor = val;
        end
    end
end

