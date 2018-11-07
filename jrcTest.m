import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin

tr = TestRunner.withTextOutput;
tr.addPlugin(CodeCoveragePlugin.forFolder(jrclust.utils.basedir))