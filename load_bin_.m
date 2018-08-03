%--------------------------------------------------------------------------
function mnWav = load_bin_(vcFile, dataType, dimm, header)
    % mnWav = load_bin_(vcFile, dimm, dataType)
    % mnWav = load_bin_(fid, dimm, dataType)
    % header: header bytes

    if nargin<2, dataType = []; end
    if nargin<3, dimm = []; end
    if nargin<4, header = 0; end
    if isempty(dataType), dataType = 'int16'; end
    mnWav = [];

    if ischar(vcFile)
        fid = [];
        if ~fileExists(vcFile)
            fprintf(2, 'File does not exist: %s\n', vcFile);
            return;
        end
        fid = fopen(vcFile, 'r');
        if header>0, fseek(fid, header, 'bof'); end
        if isempty(dimm) % read all
            S_file = dir(vcFile);
            if numel(S_file)~=1, return; end % there must be one file
            nData = floor((S_file(1).bytes - header) / bytesPerSample_(dataType));
            dimm = [nData, 1]; %return column
        end
    else % fid directly passed
        fid = vcFile;
        if isempty(dimm), dimm = inf; end
    end
    try
        t1 = tic;
        mnWav = fread_(fid, dimm, dataType);
        if ischar(vcFile)
            fclose(fid);
            fprintf('Loading %s took %0.1fs\n', vcFile, toc(t1));
        end
    catch
        disperr_();
    end
end %func
