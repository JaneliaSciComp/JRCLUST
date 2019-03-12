function fileData = readHeader(rawPath)
    %READHEADER Read in and store header data
    % adapted from scripts available here: http://www.intantech.com/RHD2000_evaluation_system.html
    d = dir(rawPath);
    if isempty(d)
        warning('Could not find file %s', rawPath);
        fileData = [];
        return;
    end

    fSizeBytes = d.bytes;
    fid = fopen(rawPath, 'r');

    % Check 'magic number' at beginning of file to make sure this is an Intan
    % Technologies RHD2000 data file.
    if fread(fid, 1, 'uint32') ~= hex2dec('c6912702')
        error('Unrecognized file type.');
    end

    % Read version number.
    fileData.versionMajor = fread(fid, 1, 'int16');
    fileData.versionMinor = fread(fid, 1, 'int16');

    if (fileData.versionMajor == 1)
        fileData.nSamplesBlock = 60;
    else
        fileData.nSamplesBlock = 128;
    end

    % Next is a block of data containing global sampling rate and amplifier frequency parameters.
    fileData.sampleRate = double(fread(fid, 1, 'single'));
    fileData.dspEnabled = fread(fid, 1, 'int16'); % DSP offset removal high-pass filter was disabled if 0, enabled if 1

    % The RHD2000 chips are not always capable of achieving the precise cutoff frequencies specified by the user,
    % so both the values requested in the GUI and the actual values realized on the chip are saved
    fileData.dspCutoffFreqActual = fread(fid, 1, 'single'); % DSP offset removal high-pass filter cutoff frequency (units: Hz)
    fileData.lowerBandwidthActual = fread(fid, 1, 'single'); % Amplifier analog high-pass filter cutoff frequency (units: Hz)
    fileData.upperBandwidthActual = fread(fid, 1, 'single'); % Amplifier analog low-pass filter cutoff frequency (units: Hz)

    fileData.dspCutoffFreqRequested = fread(fid, 1, 'single'); % User-requested DSP offset removal filter cutoff frequency (units: Hz)
    fileData.lowerBandwidthRequested = fread(fid, 1, 'single'); % User-requested amplifier high-pass filter cutoff frequency (units: Hz)
    fileData.upperBandwidthRequested = fread(fid, 1, 'single'); % User-requested amplifier low-pass filter cutoff frequency (units: Hz)

    % The next parameter records the state of the software-implemented 50/60 Hz notch filter in the GUI during recording.
    % This notch filter is never applied to saved data, but this information may be used to re-apply the notch filter
    % to recorded data, if desired.
    fileData.notchFilterMode = fread(fid, 1, 'int16'); % Software notch filter was disabled if 0...
    fileData.notchFilterFreq = 0;
    if (fileData.notchFilterMode == 1)
        fileData.notchFilterFreq = 50; % ...enabled and set to 50 Hz if 1...
    elseif (fileData.notchFilterMode == 2)
        fileData.notchFilterFreq = 60; % ...or enabled and set to 60 Hz if 2
    end

    % Next are two floating-point numbers indicating the latest user-requested electrode impedance test frequency
    % and the impedance test frequency actually realized on the RHD2000 chip
    fileData.impedanceTestFreqRequested = fread(fid, 1, 'single'); % Electrode impedance test frequency last requested by user (units: Hz)
    fileData.impedanceTestFreqActual = fread(fid, 1, 'single'); % Closest realizable electrode impedance test frequency (units: Hz)

    % In the “Configure” tab of the Intan GUI, there are three general-purpose text fields that may be used to enter notes
    % on particular recording sessions. The contents of these text fields are saved here.
    fileData.notes = {freadQString(fid), ... % User text from the “Note 1” field in the “Configure” tab in the GUI.
                      freadQString(fid), ... % User text from the “Note 2” field in the “Configure” tab in the GUI.
                      freadQString(fid)};    % User text from the “Note 3” field in the “Configure” tab in the GUI.

    % Versions 1.1 and later support saving on-chip temperature sensor readings.
    % The following number is saved in these versions, indicating the number of temperature sensors recorded.
    % This number is typically equal to the number of RHD2000 chips plugged into the SPI ports, or zero if the temperature recording option is disabled.
    % Note that some file formats do not support saving temperature sensor data; in those formats, this number will always be zero.
    fileData.nTempSensors = 0; % Number of temperature sensor channels saved in file. This number is set to zero if the option for saving temperature data is not selected.
    if ((fileData.versionMajor == 1 && fileData.versionMinor >= 1) || (fileData.versionMajor > 1))
        fileData.nTempSensors = fread(fid, 1, 'int16');
    end

    % Version 1.3 and later saves the “board mode”. This integer is set by hardware in the RHD2000 USB interface board or Intan Recording Controller.
    % Currently, a board mode of zero indicates that the on-board ADCs operate over a range of 0-3.3V;
    % a board mode of one indicates that the on-board ADCs operate over a range of ±5.0V.
    % The Intan Recording Controller has a board mode of 13. This indicates that the ADCs (Analog In ports) operate over a range of ±10.24V.
    fileData.boardMode = 0;
    if ((fileData.versionMajor == 1 && fileData.versionMinor >= 3) || (fileData.versionMajor > 1))
        fileData.boardMode = fread(fid, 1, 'int16');
    end

    % Version 2.0 and later (supported only by the 512- or 1024-channel Intan Recording Controller) saves name of the channel used for digital re-referencing.
    % The waveform from this channel may be added to other amplifier channels to undo the effects of digital re-referencing, if desired.
    % If hardware referencing was selected, this string is set to “n/a”.
    if (fileData.versionMajor > 1)
        fileData.refChannelName = freadQString(fid);
    end

    % Define data structure for spike trigger settings.
    spikeTriggerStruct = struct('voltageTriggerMode', {}, ...    % 0: trigger on digital input; 1: trigger on voltage threshold 
                                'voltageThreshold', {}, ...      % Spike voltage threshold (units: microvolts)
                                'digitalTriggerChannel', {}, ... % USB board digital input channel used for spike trigger (0-15)
                                'digitalEdgePolarity', {});      % 0: trigger on digital falling edge; 1: trigger on digital rising edge

    newTriggerChannel = struct(spikeTriggerStruct);
    spikeTriggers = struct(spikeTriggerStruct);

    % Define data structure for data channels.
    channelStruct = struct('nativeChannelName', {}, ... % e.g., “Port B” or “Board Digital Inputs”
                           'customChannelName', {}, ... % e.g., “B” or “DIN”
                           'nativeOrder', {}, ...       % The original numerical order of this channel in the GUI display (e.g., the native order of amplifier channel B-013 is 13).
                           'customOrder', {}, ...       % The numerical order of this channel as it appears on the GUI, after possible reordering by the user.
                           'boardStream', {}, ...       % USB board data stream (0-7); each data stream supports up to 32 channels. Each RHD2164 chips use two data streams. Each RHD2216 chip uses an entire data stream.
                           'chipChannel', {}, ...       % RHD2000 channel number (0-31)
                           'portName', {}, ...          % aka signalGroupName
                           'portPrefix', {}, ...        % aka signalGroupPrefix
                           'portNumber', {}, ...        % aka group ID
                           'electrodeImpedanceMagnitude', {}, ... % Last measured impedance magnitude (units: Ohms)
                           'electrodeImpedancePhase', {}); % Last measured impedance phase (units: degrees)

    newChannel = struct(channelStruct);

    % Create structure arrays for each type of data channel.
    fileData.ampChannels = struct(channelStruct);
    fileData.auxInputChannels = struct(channelStruct);
    fileData.supplyVoltageChannels = struct(channelStruct);
    fileData.boardADCChannels = struct(channelStruct);
    fileData.boardDigInChannels = struct(channelStruct);
    fileData.boardDigOutChannels = struct(channelStruct);

    ampIndex = 1;
    auxInputIndex = 1;
    supplyVoltageIndex = 1;
    boardADCIndex = 1;
    boardDigInIndex = 1;
    boardDigOutIndex = 1;

    % The next number indicates the number of “signal groups” present in the data file.
    % This number is typically equal to seven: Port A, Port B, Port C, Port D, Board ADC Inputs, Board Digital Inputs, and Board Digital Outputs.
    % If a 1024-channel Intan Recording Controller was used, this number will be equal to 11: Ports A – H, Board ADC Inputs, Board Digital Inputs, and Board Digital Outputs.
    nSignalGroups = fread(fid, 1, 'int16');

    for iGroup = 1:nSignalGroups
        % For each signal group, the following “signal group header” is saved, along with a description of each channel in the signal group.
        signalGroupName = freadQString(fid);
        signalGroupPrefix = freadQString(fid);
        signalGroupEnabled = fread(fid, 1, 'int16'); % disabled if 0, enabled if 1
        nChannelsGroup = fread(fid, 1, 'int16'); % Total number of channels in signal group
        nAmpChannelsGroup = fread(fid, 1, 'int16'); % Of the total number of channels in the signal group, the number that are amplifier channels

        % Immediately following a signal group (before the remaining signal group headers) is a list of channel descriptions.
        % If a signal group is enabled and has more than zero channels, then for each channel the following information is saved.
        if (nChannelsGroup > 0 && signalGroupEnabled > 0)
            newChannel(1).portName = signalGroupName;
            newChannel(1).portPrefix = signalGroupPrefix;
            newChannel(1).portNumber = iGroup;

            for iChannel = 1:nChannelsGroup
                newChannel(1).nativeChannelName = freadQString(fid);
                newChannel(1).customChannelName = freadQString(fid);
                newChannel(1).nativeOrder = fread(fid, 1, 'int16');
                newChannel(1).customOrder = fread(fid, 1, 'int16');

                signalType = fread(fid, 1, 'int16'); % see below

                channelEnabled = fread(fid, 1, 'int16'); % 0: channel disabled; 1: channel enabled

                newChannel(1).chipChannel = fread(fid, 1, 'int16');
                newChannel(1).boardStream = fread(fid, 1, 'int16');

                % Even non-amplifier channels will contain fields for Spike Scope trigger parameters and electrode impedance data,
                % but these fields will contain default values that may be ignored
                newTriggerChannel(1).voltageTriggerMode = fread(fid, 1, 'int16');
                newTriggerChannel(1).voltageThreshold = fread(fid, 1, 'int16');
                newTriggerChannel(1).digitalTriggerChannel = fread(fid, 1, 'int16');
                newTriggerChannel(1).digitalEdgePolarity = fread(fid, 1, 'int16');
                newChannel(1).electrodeImpedanceMagnitude = fread(fid, 1, 'single');
                newChannel(1).electrodeImpedancePhase = fread(fid, 1, 'single');

                if (channelEnabled)
                    switch (signalType)
                        case 0 % RHD2000 amplifier channel
                            fileData.ampChannels(ampIndex) = newChannel;
                            spikeTriggers(ampIndex) = newTriggerChannel;
                            ampIndex = ampIndex + 1;
                        case 1 % RHD2000 auxiliary input channel
                            fileData.auxInputChannels(auxInputIndex) = newChannel;
                            auxInputIndex = auxInputIndex + 1;
                        case 2 % RHD2000 supply voltage channel
                            fileData.supplyVoltageChannels(supplyVoltageIndex) = newChannel;
                            supplyVoltageIndex = supplyVoltageIndex + 1;
                        case 3 % USB board ADC input channel
                            fileData.boardADCChannels(boardADCIndex) = newChannel;
                            boardADCIndex = boardADCIndex + 1;
                        case 4 % USB board digital input channel
                            fileData.boardDigInChannels(boardDigInIndex) = newChannel;
                            boardDigInIndex = boardDigInIndex + 1;
                        case 5 % USB board digital output channel
                            fileData.boardDigOutChannels(boardDigOutIndex) = newChannel;
                            boardDigOutIndex = boardDigOutIndex + 1;
                        otherwise
                            error('Unknown channel type');
                    end
                end
            end
        end
    end

    fileData.headerOffset = ftell(fid);

    % Summarize contents of data file.
    fileData.nChans = ampIndex - 1;
    fileData.nAuxInputChannels = auxInputIndex - 1;
    fileData.nSupplyVoltageChannels = supplyVoltageIndex - 1;
    fileData.nBoardADCChannels = boardADCIndex - 1;
    fileData.nBoardDigInChannels = boardDigInIndex - 1;
    fileData.nBoardDigOutChannels = boardDigOutIndex - 1;

    % Determine how many samples the data file contains. Each data block contains N amplifier samples,
    % where N = 60 for the 256-channel RHD2000 Evaluation System and N = 128 for the 512-channel and 1024-channel Intan Recording Controllers. 

    % Sequential integers (e.g., 0, 1, 2, 3...) with zero marking the beginning of a recording or a trigger point.
    % Time indices can be negative to denote pre-trigger times.
    % Divide by the amplifier sampling rate (in Samples/s) to get a time vector with units of seconds.
    fileData.nBytesBlock = fileData.nSamplesBlock * 4;  % timestamp data (int32)

    % For each enabled RHD2000 amplifier channel, N ADC samples (uint16):
    fileData.nBytesBlock = fileData.nBytesBlock + fileData.nSamplesBlock * 2 * fileData.nChans;

    % For each enabled RHD2000 auxiliary input channel, N/4 ADC samples (uint16).
    % (The Intan GUI software samples the RHD2000 auxiliary inputs at one-fourth the sampling rate of the amplifiers.)
    fileData.nBytesBlock = fileData.nBytesBlock + 15 * 2 * fileData.nAuxInputChannels;

    % For each enabled RHD2000 supply voltage channel, one ADC sample (uint16).
    % (The Intan GUI software samples the RHD2000 supply voltage once per data block.)
    fileData.nBytesBlock = fileData.nBytesBlock + 1 * 2 * fileData.nSupplyVoltageChannels;

    % For each enabled RHD2000 temperature sensor channel, one ADC sample (int16).
    % (The Intan GUI software samples the RHD2000 temperature sensor once per data block.)
    % The GUI calculates a running average of temperature sensor readings over a window of approximately 100 ms to improve sensor accurate.
    if (fileData.nTempSensors > 0)
       fileData.nBytesBlock = fileData.nBytesBlock + 1 * 2 * fileData.nTempSensors; 
    end

    % For each USB interface board ADC channel, N samples (uint16):
    fileData.nBytesBlock = fileData.nBytesBlock + fileData.nSamplesBlock * 2 * fileData.nBoardADCChannels;

    % If ANY USB interface board digital inputs are enabled, unsigned 16-bit integers record N samples
    % from ALL digital inputs 0-15. If no digital inputs are enabled, these samples are not recorded.
    if (fileData.nBoardDigInChannels > 0)
        fileData.nBytesBlock = fileData.nBytesBlock + fileData.nSamplesBlock * 2;
    end

    % Board digital outputs are sampled at same rate as amplifiers
    if (fileData.nBoardDigOutChannels > 0)
        fileData.nBytesBlock = fileData.nBytesBlock + fileData.nSamplesBlock * 2;
    end

    bytesRemaining = fSizeBytes - fileData.headerOffset;
    fileData.nDataBlocks = bytesRemaining / fileData.nBytesBlock;
    fileData.nAmpSamples = fileData.nSamplesBlock * fileData.nDataBlocks;

    fileData.dataPresent = (bytesRemaining > 0);
    fileData.nAuxInputSamples = (fileData.nSamplesBlock/4) * fileData.nDataBlocks;
    fileData.nSupplyVoltageSamples = fileData.nDataBlocks;
    fileData.nBoardADCSamples = fileData.nSamplesBlock * fileData.nDataBlocks;
    fileData.nBoardDigInSamples = fileData.nSamplesBlock * fileData.nDataBlocks;
    fileData.nBoardDigOutSamples = fileData.nSamplesBlock * fileData.nDataBlocks;

    fileData.recordingTime = fileData.nAmpSamples / fileData.sampleRate;
end

%% LOCAL FUNCTIONS
function str = freadQString(fid)
    % str = freadQString(fid)
    %
    % Read Qt style QString.  The first 32-bit unsigned number indicates
    % the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
    % the string is null.

    str = '';
    length = fread(fid, 1, 'uint32');
    if length == hex2num('ffffffff')
        return;
    end
    % convert length from bytes to 16-bit Unicode words
    length = length / 2;

    for i=1:length
        str(i) = fread(fid, 1, 'uint16');
    end
end

function out = notchFilter(in, fSample, fNotch, Bandwidth)

    % out = notchFilter(in, fSample, fNotch, Bandwidth)
    %
    % Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
    % fSample = sample rate of data (in Hz or Samples/sec)
    % fNotch = filter notch frequency (in Hz)
    % Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
    %   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
    %   poor time-domain properties with an extended ringing response to
    %   transient disturbances.
    %
    % Example:  If neural data was sampled at 30 kSamples/sec
    % and you wish to implement a 60 Hz notch filter:
    %
    % out = notchFilter(in, 30000, 60, 10);

    tstep = 1/fSample;
    Fc = fNotch*tstep;

    L = length(in);

    % Calculate IIR filter parameters
    d = exp(-2*pi*(Bandwidth/2)*tstep);
    b = (1 + d*d)*cos(2*pi*Fc);
    a0 = 1;
    a1 = -b;
    a2 = d*d;
    a = (1 + d*d)/2;
    b0 = 1;
    b1 = -2*cos(2*pi*Fc);
    b2 = 1;

    out = zeros(size(in));
    out(1) = in(1);  
    out(2) = in(2);
    % (If filtering a continuous data stream, change out(1) and out(2) to the
    %  previous final two values of out.)

    % Run filter
    for i=3:L
        out(i) = (a*b2*in(i-2) + a*b1*in(i-1) + a*b0*in(i) - a2*out(i-2) - a1*out(i-1))/a0;
    end

end
