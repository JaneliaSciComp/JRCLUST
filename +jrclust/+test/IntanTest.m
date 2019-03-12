classdef IntanTest < matlab.unittest.TestCase
    %INTANTEST Test Intan recordings
    properties
        hJRC;
        intanRes;
    end

    properties (Dependent)
        hCfg;
    end
    
    %% SETUP
    methods (TestClassSetup)
        function setupProps(obj)
            hCfg_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'intan', 'test.prm'));
            obj.hJRC = JRC(hCfg_);
        end

        function getIntanRes(obj)
            obj.intanRes = cell(numel(obj.hCfg.rawRecordings), 1);
            for i = 1:numel(obj.intanRes)
                filename = obj.hCfg.rawRecordings{i};
                obj.intanRes{i} = jrclust.test.read_Intan_RHD2000_file(filename);
            end
        end
    end

    %% TEST METHODS
    methods (Test)
        function testSizes(obj)
            %TESTSIZES Assert all ampData sizes are the same
            for i = 1:numel(obj.hCfg.rawRecordings)
                hRec = jrclust.detect.newRecording(obj.hCfg.rawRecordings{i}, obj.hCfg);
                obj.assertEqual(size(obj.intanRes{i}, 1), hRec.nChans);
                obj.assertEqual(size(obj.intanRes{i}, 2), hRec.nSamples);
            end
        end

        function testROIs(obj)
            for i = 1:numel(obj.hCfg.rawRecordings)
                hRec = jrclust.detect.newRecording(obj.hCfg.rawRecordings{i}, obj.hCfg);
                rows = randsample(hRec.nChans, floor(hRec.nChans/2));
                cols = randsample(hRec.nSamples, min(hRec.nSamples, 10191));

                iRes = single(obj.intanRes{i});
                obj.assertEqual(iRes(rows, cols), hRec.readRawROI(rows, cols));
            end
        end
    end

    %% GETTERS/SETTERS
    methods
        function hCfg = get.hCfg(obj)
            hCfg = obj.hJRC.hCfg;
        end
        function set.hCfg(obj, hCfg)
            obj.hJRC.hCfg = hCfg;
        end
    end
end

