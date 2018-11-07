function abstr = about()
    %ABOUT Get an about string
    md = jrclust.utils.info();
    abstr = strjoin({jrclust.utils.version();
                     sprintf('  Last updated on %s', md.changeDate);
                     '  Created by James Jun';
                     '  With contributions from:';
                     '    Mike Economo';
                     '    Hidehiko Inagaki';
                     '    Tim Wang';
                     '  Maintained by:';
                     '    Alan Liddell, Vidrio Technologies';
                     '';
                     'Hardware Requirements';
                     '  32GB ram (or 1/4 of recording size)';
                     '  Nvidia GPU (Compute Capability 3.5+: Kepler, Maxwell or Pascal)';
                     '';
                     'Software Requirements';
                     '  Matlab (R2014b or higher) with the following toolboxes:';
                     '    Parallel Processing';
                     '    Image processing';
                     '    Signal Processing';
                     '    Statistics and Machine Learning';
                     '  CUDA version supported by Parallel Processing toolbox:';
                     '    CUDA 8.0 (R2017a,b)';
                     '    CUDA 7.5 (R2016a,b)';
                     '    CUDA 7.0 (R2015b)';
                     '    CUDA 6.5 (R2015a)';
                     '    CUDA 6.0 (R2014b)';
                     '    All CUDA versions may be found at: https://developer.nvidia.com/cuda-toolkit-archive';
                     '  Latest Nvidia GPU driver (link below)';
                     '    http://www.nvidia.com/Download/index.aspx';
                     '  Visual Studio 2015 Community (link below)';
                     '    https://docs.microsoft.com/en-us/visualstudio/install/install-visual-studio-2015?view=vs-2015'
                     }, '\n');
end