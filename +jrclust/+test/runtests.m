function res = runtests()
    close all; clear;
    res = struct();

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.CodeCoveragePlugin

    tr = TestRunner.withTextOutput;
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust'));
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.detect'));
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.features'));
    tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.sort'));
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

    detectMultiSuite = TestSuite.fromClass(?jrclust.test.DetectMultiTest);
    detectMultiRes = tr.run(detectMultiSuite);

    res = jrclust.utils.mergeStructs(res, detectMultiRes);

    %% TEST SORTING
    sortSuite = TestSuite.fromClass(?jrclust.test.SortTest);
    sortRes = tr.run(sortSuite);

    res = jrclust.utils.mergeStructs(res, sortRes);

    %% TEST CURATION
    curateSuite = TestSuite.fromClass(?jrclust.test.ManualTest);
    curateRes = tr.run(curateSuite);

    res = jrclust.utils.mergeStructs(res, curateRes);

    %% TEST INTAN
    intanSuite = TestSuite.fromClass(?jrclust.test.IntanTest);
    intanRes = tr.run(intanSuite);

    res = jrclust.utils.mergeStructs(res, intanRes);
end