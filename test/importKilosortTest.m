function tests = importKilosortTest
%IMPORTKILOSORTTEST tests for import_ksort_
tests = functiontests(localfunctions);
end % function

% SETUP
function setupOnce(testCase)  % do not change function name
dataDir = getenv('JRCTESTDATA');
fatalAssertNotEmpty(testCase, dataDir);
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
function testTimesTemplatesSites(testCase)
dataDir = getenv('JRCTESTDATA');
testDir = fullfile(dataDir, 'importKilosort', 'eMouse');
% load up test results
S0 = load(fullfile(testDir, 'test_jrc.mat'));

dataSources = fullfile(fileparts(dataDir), 'data_sources');
rezFile = fullfile(dataSources, 'eMouse', 'rez.mat');

% load up the rez file
load(rezFile, 'rez');
% lifted  from rezToPhy
spikeTimes = uint64(rez.st3(:,1));
spikeTemplates = uint32(rez.st3(:,2));

assertTrue(testCase, all(spikeTimes == S0.viTime_spk));
assertTrue(testCase, all(spikeTemplates == S0.S_clu.viTemplate_spk));
assertTrue(testCase, all(spikeTemplates == S0.S_clu.viClu_auto));
end

% TEARDOWN
function teardownOnce(testCase)
dataDir = getenv('JRCTESTDATA');
testDir = fullfile(dataDir, 'importKilosort', 'eMouse');
prmFile = fullfile(testDir, 'test.prm');

jrc('clear', prmFile);

% clean up files
delete(fullfile(testDir, 'test.prm'));
delete(fullfile(testDir, 'test_full.prm'));
delete(fullfile(testDir, 'test_spkfet.jrc'));
delete(fullfile(testDir, 'test_spkraw.jrc'));
delete(fullfile(testDir, 'test_spkwav.jrc'));
delete(fullfile(testDir, 'test-probe.mat'));
delete(fullfile(testDir, 'test_jrc.mat'));
end