classdef UtilsTest < matlab.unittest.TestCase
    %UTILSTEST Test that various utilities run as expected

    properties
        testDir;
        testFile;
        testText;
    end

    methods (TestClassSetup)
        function createTextFile(testCase)
            testCase.testDir = tempname();
            mkdir(testCase.testDir);

            testCase.testFile = fullfile(testCase.testDir, 'test.txt');
        end

        function createTestFiles(testCase)
            % make another file called test.txt in a different directory
            mkdir(testCase.testDir, 'dupedir');
            fid = fopen(fullfile(testCase.testDir, 'dupedir', 'test.txt'), 'w');
            fwrite(fid, 'aaa');
            fclose(fid);

            % make an empty directory
            mkdir(testCase.testDir, 'emptydir');
        end

        function lipsum(testCase)
            % https://tunaipsum.com/
            testCase.testText = {'Pacific herring, "driftfish cownose ray." Whalefish giant wels toadfish; eel pencilsmelt barred danio pilot fish elephant fish yellowbelly tail catfish', ...
                                 'regal whiptail catfish whitetip reef shark. Speckled trout North American darter loach; killifish mako shark. Beardfish California halibut wasp fish', ...
                                 'orbicular batfish lightfish loach catfish mosquitofish searobin rock bass cutthroat trout. Razorfish croaker, slimy mackerel topminnow jack', ...
                                 'scabbard fish Devario flier zebrafish splitfin, ground shark pickerel, bleak goosefish eel cod. Bangus cichlid snipefish monkfish pelican eel', ...
                                 'featherback paperbone scat deepwater cardinalfish. Butterflyfish coffinfish rockling redtooth triggerfish wallago grunion temperate bass. Pearlfish', ...
                                 'Pacific hake redtooth triggerfish mail-cheeked fish yellow moray herring smelt sand knifefish European minnow halosaur northern squawfish', ...
                                 'Raccoon butterfly fish. Yellowtail Ragfish zebra pleco dogteeth tetra round stingray cherry salmon Australian grayling. Crestfish: sergeant major', ...
                                 'ghost flathead, sleeper shark, man-of-war fish oceanic flyingfish velvet-belly shark, righteye flounder dogfish shark pearlfish pink salmon. Pilot fish,', ...
                                 'loach catfish common carp.'};
            fid = fopen(testCase.testFile, 'w');
            for i = 1:numel(testCase.testText)
                fprintf(fid, '%s\n', testCase.testText{i});
            end
            fclose(fid);
        end
    end

    methods (TestClassTeardown)
        function clearTestFiles(testCase)
            try
                rmdir(testCase.testDir, 's');
            catch % temp dir, it's okay if we fail
            end
        end
    end

    methods (Test)
        function about(testCase)
            aboutstr = jrclust.utils.about();
            testCase.assertSubstring(aboutstr, 'NVIDIA GPU, compute capability 3.5 or later');
        end

        function absPath(testCase)
            cd(testCase.testDir); % contains test.txt

            % give absolute path to file
            p = jrclust.utils.absPath(testCase.testFile);
            testCase.assertEqual(p, testCase.testFile);

            % give relative path to file in current directory, no hint
            p = jrclust.utils.absPath('test.txt');
            testCase.assertEqual(p, testCase.testFile);

            % give relative path to file in current directory, no hint
            p = jrclust.utils.absPath('test.txt', testCase.testDir);
            testCase.assertEqual(p, testCase.testFile);

            cd(fullfile(testCase.testDir, 'emptydir')); % empty directory

            % give relative path to file in another directory, no hint
            p = jrclust.utils.absPath('test.txt');
            testCase.assertEmpty(p);

            % give relative path to file in another directory, with hint
            p = jrclust.utils.absPath('test.txt', testCase.testDir);
            testCase.assertEqual(p, testCase.testFile);

            cd(testCase.testDir);

            % give relative path to file with same name in different
            % directory, with hint
            p = jrclust.utils.absPath('test.txt', 'dupedir');
            testCase.assertEqual(p, fullfile(testCase.testDir, 'dupedir', 'test.txt'));

            p = jrclust.utils.absPath('dupedir/test.txt');
            testCase.assertEqual(p, fullfile(testCase.testDir, 'dupedir', 'test.txt'));

            p = jrclust.utils.absPath('test.txt', fullfile(testCase.testDir, 'dupedir'));
            testCase.assertEqual(p, fullfile(testCase.testDir, 'dupedir', 'test.txt'));
        end

        function basedir(testCase)
            bd = jrclust.utils.basedir();
            testCase.assertEqual(bd, fileparts(fileparts(mfilename('fullpath'))));
        end

        function readLines(testCase)
            lines = jrclust.utils.readLines(testCase.testFile);
            testCase.fatalAssertEqual(numel(lines), numel(testCase.testText));
            for i = 1:numel(testCase.testText)
                testCase.assertEqual(lines{i}, testCase.testText{i});
            end
        end
    end
end

