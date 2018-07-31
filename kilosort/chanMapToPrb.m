%--------------------------------------------------------------------------
function viSiteZero = chanMapToPrb(chanMap, prbFilename)
    if ischar(chanMap)
        chanMap = load(chanMap);
    end
    
    fid = fopen(prbFilename, 'w');
    
    viSiteZero = find(~chanMap.connected)';
    
    channels = chanMap.chanMap;
    fprintf(fid, 'channels = [%s];\n', strjoin(strsplit(num2str(channels(:)'))));
    fprintf(fid, 'xcoords = [%s];\n', strjoin(strsplit(num2str(chanMap.xcoords(:)'))));
    fprintf(fid, 'ycoords = [%s];\n', strjoin(strsplit(num2str(chanMap.ycoords(:)'))));
    fprintf(fid, 'geometry = [xcoords'' ycoords''];\n');
    fprintf(fid, 'pad = [12 12];\n'); % verify
    
    if isfield(chanMap, 'kcoords')
        fprintf(fid, 'shank = [%s];\n', strjoin(strsplit(num2str(chanMap.kcoords(:)'))));
    elseif isfield(chanMap, 'shankInd')
        fprintf(fid, 'shank = [%s];\n', strjoin(strsplit(num2str(chanMap.shankInd(:)'))));
    end
    
    fclose(fid);    
end %func