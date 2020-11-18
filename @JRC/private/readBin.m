function binData = readBin(filename, binShape, dataType)
%LOADBIN Load traces/features from binary file
if exist(filename, 'file')
    fid = fopen(filename, 'r');
    binData = fread(fid, Inf, dataType);
    fclose(fid);
    binData = reshape(binData, binShape);
else
    binData = [];
end
end %fun