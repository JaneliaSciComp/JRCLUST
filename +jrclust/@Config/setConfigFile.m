function setConfigFile(obj, configFile, reloadParams)
    %SETCONFIGFILE Don't use this. Seriously.
     if nargin < 2
         return;
     end
     if nargin < 3
         reloadParams = 1;
     end

     configFile_ = jrclust.utils.absPath(configFile);
     if isempty(configFile_)
         error('Could not find %s', configFile);
     end

     obj.configFile = configFile_;

     if reloadParams
         obj.loadParams(obj.configFile);
     end
end