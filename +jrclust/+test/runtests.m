function res = runtests()
close all; clear;

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin

tr = TestRunner.withTextOutput;
tr.addPlugin(CodeCoveragePlugin.forFolder([jrclust.utils.basedir '/@JRC'], ...
    'IncludingSubfolders', true));
% tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust'));
tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.interfaces'));
% tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.detect'));
% tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.features'));
tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.sort'));
% tr.addPlugin(CodeCoveragePlugin.forPackage('jrclust.curate'));

suites = [];

%% TEST DELETE
suites = [suites TestSuite.fromClass(?jrclust.test.DensityPeakClustering.DeleteTest)];

%% TEST UNDELETE
suites = [suites TestSuite.fromClass(?jrclust.test.DensityPeakClustering.UndeleteTest)];

%% TEST MERGE
suites = [suites TestSuite.fromClass(?jrclust.test.DensityPeakClustering.MergeTest)];

%% TEST SPLIT
suites = [suites TestSuite.fromClass(?jrclust.test.DensityPeakClustering.SplitTest)];

%% TEST REVERT
suites = [suites TestSuite.fromClass(?jrclust.test.DensityPeakClustering.RevertTest)];

%% TEST CONVERT HISTORY
suites = [suites TestSuite.fromClass(?jrclust.test.DensityPeakClustering.ConvertHistoryTest)];

%% RUN TESTS, SUMMARIZE
res = tr.run(suites);

if ~any([res.Failed])
    disp('all tests passed! :^)');
end
end % func