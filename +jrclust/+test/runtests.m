function res = runTests()
    close all; clear;
    res = struct();

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.CodeCoveragePlugin

    tr = TestRunner.withTextOutput;
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust'));
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.controllers.detect'));

    %% TEST BOOTSTRAP
%     bootstrapSuite = TestSuite.fromClass(?jrclust.test.BootstrapTest);
%     bootstrapRes = tr.run(bootstrapSuite);
% 
%     res = jrclust.utils.mergeStructs(res, bootstrapRes);

    %% TEST DETECTION
    detectSuite = TestSuite.fromClass(?jrclust.test.DetectTest);
    detectRes = tr.run(detectSuite);

    res = jrclust.utils.mergeStructs(res, detectRes);
end