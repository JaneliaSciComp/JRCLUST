function tests = manualKilosortTest
%MANUALKILOSORTTEST tests for manually curating an imported KiloSort
%session
tests = functiontests(localfunctions);
end

% SETUP
function setupOnce(testCase)
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
% fire up the gui
jrc('manual', prmFile);
end

% TESTS
function testDeleteCluster(testCase)
S0 = get0_();
S_clu_before = S0.S_clu;
simsBefore = S_clu_before.mrSim_clu;
tempsBefore = S_clu_before.cviTemplate_clu;
spikesBefore = S_clu_before.cviSpk_clu;

% delete cluster 1
S0.iCluCopy = 1;
S0 = ui_delete_(S0);
S_clu_after = S0.S_clu;
simsAfter = S_clu_after.mrSim_clu;
tempsAfter = S_clu_after.cviTemplate_clu;
spikesAfter = S_clu_after.cviSpk_clu;

% test sim scores
% sim scores in unaffected clusters are unchanged
fatalAssertTrue(testCase, all(size(simsAfter) == size(simsBefore) - 1));
assertEqual(testCase, norm(simsBefore(2:end, 2:end) - simsAfter), 0);
        
% test templates
% templates in unaffected clusters are unchanged
fatalAssertEqual(testCase, numel(tempsAfter), numel(tempsBefore) - 1);
assertTrue(testCase, all(cellfun(@(x, y) all(x == y), ...
           tempsBefore(2:end), tempsAfter)));
       
% test spikes
% spikes in unaffected clusters are unchanged
fatalAssertEqual(testCase, numel(spikesAfter), numel(spikesBefore) - 1);
assertTrue(testCase, all(cellfun(@(x, y) all (x == y), ...
           spikesBefore(2:end), spikesAfter)));
end

function testSplitCluster(testCase)
S0 = get0_();
S_clu_before = S0.S_clu;
simsBefore = S_clu_before.mrSim_clu;
tempsBefore = S_clu_before.cviTemplate_clu;
spikesBefore = S_clu_before.cviSpk_clu;

% split cluster 1 clean in half
clusterSpikes = spikesBefore{1};
vlIn = true(size(clusterSpikes));
vlIn(floor(numel(clusterSpikes) / 2):end) = 0;

S_clu_after = split_clu_(1, vlIn);
simsAfter = S_clu_after.mrSim_clu;
tempsAfter = S_clu_after.cviTemplate_clu;
spikesAfter = S_clu_after.cviSpk_clu;

% test sim scores
% sim scores in unaffected clusters are unchanged
fatalAssertTrue(testCase, ...
                all(size(simsAfter) == size(simsBefore) + 1));
assertEqual(testCase, ...
            norm(simsBefore(2:end, 2:end) - simsAfter(3:end, 3:end)), 0);
        
% test templates
% templates in unaffected clusters are unchanged
fatalAssertEqual(testCase, numel(tempsAfter), numel(tempsBefore) + 1);
assertTrue(testCase, all(cellfun(@(x, y) all(x == y), ...
           tempsBefore(2:end), tempsAfter(3:end))));
% templates in initial cluster are distributed over new clusters
fatalAssertEqual(testCase, numel(union(tempsAfter{1}, tempsAfter{2})), ...
                 numel(tempsBefore{1}));
assertTrue(testCase, all(union(tempsAfter{1}, tempsAfter{2}) == tempsBefore{1}));

       
% test spikes
% spikes in unaffected clusters are unchanged
fatalAssertEqual(testCase, numel(spikesAfter), numel(spikesBefore) + 1);
assertTrue(testCase, all(cellfun(@(x, y) all (x == y), ...
           spikesBefore(2:end), spikesAfter(3:end))));
% spikes in initial cluster are distributed over new clusters
fatalAssertEqual(testCase, numel([spikesAfter{1}; spikesAfter{2}]), ...
                 numel(spikesBefore{1}));
assertTrue(testCase, all(sort([spikesAfter{1}; spikesAfter{2}]) == spikesBefore{1}));
end

function testMergeClusters(testCase)
S0 = get0_();
S_clu_before = S0.S_clu;
simsBefore = S_clu_before.mrSim_clu;
tempsBefore = S_clu_before.cviTemplate_clu;
spikesBefore = S_clu_before.cviSpk_clu;

% merge clusters 1 and 2
S0.iCluCopy = 1;
S0.iCluPaste = 2;

S0 = ui_merge_(S0);
S_clu_after = S0.S_clu;
simsAfter = S_clu_after.mrSim_clu;
tempsAfter = S_clu_after.cviTemplate_clu;
spikesAfter = S_clu_after.cviSpk_clu;

% test sim scores
% sim scores in unaffected clusters are unchanged
fatalAssertTrue(testCase, ...
                all(size(simsAfter) == size(simsBefore) - 1));
assertEqual(testCase, ...
            norm(simsAfter(2:end, 2:end) - simsBefore(3:end, 3:end)), 0);
        
% test templates
% templates in unaffected clusters are unchanged
fatalAssertEqual(testCase, numel(tempsAfter), numel(tempsBefore) - 1);
assertTrue(testCase, all(cellfun(@(x, y) all(x == y), ...
           tempsBefore(3:end), tempsAfter(2:end))));
% templates in merging clusters are now in merged cluster
fatalAssertEqual(testCase, numel(union(tempsBefore{1}, tempsBefore{2})), ...
                 numel(tempsAfter{1}));
assertTrue(testCase, all(union(tempsBefore{1}, tempsBefore{2}) == tempsAfter{1}));

       
% test spikes
% spikes in unaffected clusters are unchanged
fatalAssertEqual(testCase, numel(spikesAfter), numel(spikesBefore) - 1);
assertTrue(testCase, all(cellfun(@(x, y) all (x == y), ...
           spikesBefore(3:end), spikesAfter(2:end))));
% spikes in initial cluster are distributed over new clusters
fatalAssertEqual(testCase, numel(union(spikesBefore{1}, spikesBefore{2})), ...
                 numel(spikesAfter{1}));
assertTrue(testCase, all(sort([spikesBefore{1}; spikesBefore{2}]) == spikesAfter{1}));
end

% TEARDOWN
function teardownOnce(testCase)
% close the gui
exit_manual_(get_fig_cache_('FigWav'));

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
