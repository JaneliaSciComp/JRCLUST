function res = runtests()
    close all; clear;
    res = struct();

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.CodeCoveragePlugin

    tr = TestRunner.withTextOutput;
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust'));
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.detect'));
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.curate'));

    %% TEST BOOTSTRAP
%     bootstrapSuite = TestSuite.fromClass(?jrclust.test.BootstrapTest);
%     bootstrapRes = tr.run(bootstrapSuite);
% 
%     res = jrclust.utils.mergeStructs(res, bootstrapRes);

    %% TEST DETECTION
    detectSuite = TestSuite.fromClass(?jrclust.test.DetectTest);
    detectRes = tr.run(detectSuite);
    
    res = jrclust.utils.mergeStructs(res, detectRes);

    %% TEST CURATION
    curateSuite = TestSuite.fromClass(?jrclust.test.ManualTest);
    curateRes = tr.run(curateSuite);

    res = jrclust.utils.mergeStructs(res, curateRes);
end