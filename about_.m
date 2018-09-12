%--------------------------------------------------------------------------
% function about_()
% 7/24/17 JJJ: Updated requirements and contact info
function csAbout = about_(varargin)
    [vcVer, vcDate] = jrc_version_();
    csAbout = { ...
    '';
    sprintf('JRCLUST %s (jrc.m)', vcVer);
    sprintf('  Last updated on %s', vcDate);
    '  Created by James Jun (jamesjun@gmail.com)';
    '  With contributions from:';
    '    Mike Economo';
    '    Hidehiko Inagaki';
    '    Tim Wang';
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
    };
    if nargout==0, disp_cs_(csAbout); end
end
