classdef UtilsTest < matlab.unittest.TestCase
    %UTILSTEST Test that various utilities run as expected
    
    properties
        testFile;
        testText;
    end
    
    methods(TestMethodSetup)
        function createTextFile(testCase)
            testCase.testFile = [tempname '.txt'];
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
 
    methods(TestMethodTeardown)
        function baleetTextFile(testCase)
            delete(testCase.testFile);
        end
    end

    
    methods (Test)
        function about(testCase)
            aboutstr = jrclust.utils.about();
            testCase.assertSubstring(aboutstr, 'Nvidia GPU (Compute Capability 3.5+: Kepler, Maxwell or Pascal)');
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

