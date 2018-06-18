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
    '';
    'Hardware Requirements';
    '  32GB ram (or 1/4 of recording size)';
    '  NVIDIA GPU (Compute Capability 3.5+: Kepler, Maxwell or Pascal)';
    '';
    'Software Requirements';
    '  Matlab (R2014b or higher) with Toolboxes below';
    '    Parallel Processing, Image processing, Signal Processing, Statistics and Machine Learning';
    '  CUDA version supported by Matlab prallel processing toolbox (link below)';
    '    CUDA 8.0 (R2017a,b) link: https://developer.nvidia.com/cuda-downloads';
    '    CUDA 7.5 (R2016a,b) link: https://developer.nvidia.com/cuda-75-downloads-archive';
    '    CUDA 7.0 (R2015b) link: https://developer.nvidia.com/cuda-toolkit-70';
    '    CUDA 6.5 (R2015a) link: https://developer.nvidia.com/cuda-toolkit-65';
    '    CUDA 6.0 (R2014b) link: https://developer.nvidia.com/cuda-toolkit-60';
    '  Latest NVidia GPU driver (link below)';
    '    http://www.nvidia.com/Download/index.aspx';
    '  Visual Studio 2013 Express (link below)';
    '    https://www.dropbox.com/s/jrpmto1mwdp4uga/en_visual_studio_express_2013_for_windows_desktop_with_update_5_x86_web_installer_6815514.exe?dl=1'
    };
    if nargout==0, disp_cs_(csAbout); end
end
