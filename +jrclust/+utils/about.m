function abstr = about()
    %ABOUT Get an about string
    md = jrclust.utils.info();
    abstr = strjoin({sprintf('%s v%s', md.program, jrclust.utils.version());
                     sprintf('Created by %s', md.authors{1});
                     'Maintained by:';
                     sprintf('    %s', strjoin(md.maintainers, ', '));
                     '';
                     'With contributions from:';
                     sprintf('    - %s', strjoin(md.contributors, '\n    - '));
                     '';
                     'Hardware Requirements:';
                     sprintf('    - %s', strjoin(md.hardwareRequirements, '\n    - '));
                     '';
                     'Software Requirements:';
                     sprintf('    - %s', strjoin(md.softwareRequirements, '\n    - '));
                     '';
                     sprintf('Last update: %s', md.changeDate);
                     }, '\n');
end