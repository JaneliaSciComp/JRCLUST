function res = runtests(packages)
close all;

if nargin < 1
    packages = {'JRC', 'DensityPeakClustering', 'TemplateClustering', ...
                'CurateController', 'utils'};
elseif isa(packages, 'char')
    packages = {packages};
end

%% INITIALIZE TEST RUNNER
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.CodeCoveragePlugin;

testRunner = TestRunner.withTextOutput;

%% COLLECT SUITES
suites = [];
if ismember('JRC', packages)
    suites = [suites collectJRCTests()];
    testRunner.addPlugin(CodeCoveragePlugin.forFolder([jrclust.utils.basedir '/@JRC'], ...
        'IncludingSubfolders', true));
end

if ismember('DensityPeakClustering', packages)
    suites = [suites collectDensityPeakClusteringTests()];
    testRunner.addPlugin(CodeCoveragePlugin.forPackage('jrclust.interfaces'));
    testRunner.addPlugin(CodeCoveragePlugin.forPackage('jrclust.sort'));
end

if ismember('TemplateClustering', packages)
    suites = [suites collectTemplateClusteringTests()];
    testRunner.addPlugin(CodeCoveragePlugin.forPackage('jrclust.interfaces'));
    testRunner.addPlugin(CodeCoveragePlugin.forPackage('jrclust.sort'));
end

if ismember('CurateController', packages)
    suites = [suites collectCurateControllerTests()];
    testRunner.addPlugin(CodeCoveragePlugin.forPackage('jrclust.curate'));
end

if ismember('utils', packages)
    suites = [suites collectUtilsTests()];
    testRunner.addPlugin(CodeCoveragePlugin.forPackage('jrclust.utils'));
end

%% RUN TESTS, SUMMARIZE
res = testRunner.run(suites);

if ~any([res.Failed])
    disp('all tests passed! :^)');
end
end %fun

%% LOCAL FUNCTIONS
function suites = collectJRCTests()
%COLLECTJRCTESTS Collect all test suites for JRC objects.
import matlab.unittest.TestSuite;

suites = TestSuite.fromClass(?jrclust.test.JRC.ConvertHistoryTest);
end %fun

function suites = collectDensityPeakClusteringTests()
%COLLECTDENSITYPEAKCLUSTERINGTESTS Collect all test suites for
%DensityPeakClustering objects.
import matlab.unittest.TestSuite;

suites = [ ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.AutoMergeTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.DeleteTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.UndeleteTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.MergeTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.SplitTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.RevertTest) ...
];
end %fun

function suites = collectTemplateClusteringTests()
%COLLECTEMPLATECLUSTERINGTESTS Collect all test suites for
%TemplateClustering objects.
import matlab.unittest.TestSuite;

suites = [ ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.DeleteTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.UndeleteTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.MergeTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.SplitTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.RevertTest) ...
];
end %fun

function suites = collectCurateControllerTests()
%COLLECTCURATECONTROLLERTESTS Collect all test suites for CurateController
%objects.
import matlab.unittest.TestSuite;

suites = [ ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.CurateController.DeleteTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.CurateController.MergeTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.CurateController.SplitTest) ...
    TestSuite.fromClass(?jrclust.test.DensityPeakClustering.CurateController.RevertTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.CurateController.DeleteTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.CurateController.MergeTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.CurateController.SplitTest) ...
    TestSuite.fromClass(?jrclust.test.TemplateClustering.CurateController.RevertTest) ...
];
end %fun

function suites = collectUtilsTests()
import matlab.unittest.TestSuite;

suites = [ ...
    TestSuite.fromClass(?jrclust.test.utils.SaveStructTest) ...
];
end %fun