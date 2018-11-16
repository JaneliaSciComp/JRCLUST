classdef DetectTest < matlab.unittest.TestCase
    %DETECTTEST

    properties
        dRes;
        testConfig;
    end

    methods (TestClassSetup)
        function doDetect(testCase)
            testCase.testConfig = fullfile(getenv('JRCTESTDATA'), 'test-dev', 'testset.prm');
            hJRC = jrc('detect', testCase.testConfig);

            testCase.dRes = hJRC.dRes;
        end
    end

    methods (Test)
        function obj = DetectionTest(testCase)
            obj.Property1 = inputArg1 + inputArg2;
        end
    end
end

