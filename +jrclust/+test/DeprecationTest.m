classdef DeprecationTest < matlab.unittest.TestCase
    %DEPRECATIONTESTS Test that deprecated commands are indeed deprecated

    methods (Test)
        function compileKsort(testCase)
            hJRC = jrc('compile-ksort');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function doc(testCase)
            hJRC = jrc('doc');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function docEdit(testCase)
            hJRC = jrc('doc-edit');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function download(testCase)
            hJRC = jrc('download');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function gitPull(testCase)
            hJRC = jrc('git-pull');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function gui(testCase)
            hJRC = jrc('gui');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function importKilosortSort(testCase)
            hJRC = jrc('import-kilosort-sort');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function importKsortSort(testCase)
            hJRC = jrc('import-ksort-sort');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function install(testCase)
            hJRC = jrc('install');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function issue(testCase)
            hJRC = jrc('issue');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function kilosort(testCase)
            hJRC = jrc('kilosort');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function kilosortVerify(testCase)
            hJRC = jrc('kilosort-verify');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function ksort(testCase)
            hJRC = jrc('ksort');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function ksortVerify(testCase)
            hJRC = jrc('ksort-verify');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function prmSet(testCase)
            hJRC = jrc('set');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function setPrm(testCase)
            hJRC = jrc('set-prm');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function setPrm2(testCase)
            hJRC = jrc('setprm');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function update(testCase)
            hJRC = jrc('update');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function which(testCase)
            hJRC = jrc('which');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function wiki(testCase)
            hJRC = jrc('wiki');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end

        function wikiDownload(testCase)
            hJRC = jrc('wiki-download');
            testCase.assertMatches(hJRC.invocationType(), '^info$');
        end
    end
end

