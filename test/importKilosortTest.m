function tests = importKilosortTest
%IMPORTKILOSORTTEST tests for import_ksort_
if isempty(getenv('JRCTESTDATA'))
    assert(0 == 1);
else
    tests = functiontests(localfunctions);
end
end % function

% SETUP
function setupOnce(testCase)  % do not change function name
dataDir = getenv('JRCTESTDATA');
dataSources = fullfile(fileparts(dataDir), 'data_sources');
testDir = fullfile(dataDir, 'importKilosort', 'eMouse');

% write test prm file
prmFile = fullfile(testDir, 'test.prm');
fid = fopen(prmFile, 'w');
fprintf(fid, 'vcFile = ''%s'';', fullfile(dataSources, 'eMouse', 'sim_binary.dat'));
fprintf(fid, 'vcFile_rez = ''%s'';', fullfile(dataSources, 'eMouse', 'rez.mat'));
fclose(fid);

% do the import
jrc('import-ksort', prmFile);
end

% TESTS
function testFoo(testCase)
one = 1;
assertEqual(testCase, one, 1);
end

% TEARDOWN
function teardownOnce(testCase)
dataDir = getenv('JRCTESTDATA');
testDir = fullfile(dataDir, 'importKilosort', 'eMouse');

% clean up files
delete(fullfile(testDir, 'test.prm'));
delete(fullfile(testDir, 'test_spkfet.jrc'));
delete(fullfile(testDir, 'test_spkraw.jrc'));
delete(fullfile(testDir, 'test_spkwav.jrc'));
delete(fullfile(testDir, 'test-probe.mat'));
end