function [data, loadTime] = readBin(filename, dtype, dshape, offset)
    %READBIN Read a binary file into memory

    if nargin < 2 || isempty(dtype)
        dtype = 'int16';
    end
    bpt = jrclust.utils.typeBytes(dtype); % #bytes in sample type
    assert(bpt > 0, 'dtype not recognized: %s', dtype);

    if nargin < 3
        dshape = [];
    end
    assert(isnumeric(dshape) && all(dshape > 0), 'dshape not valid: %s', dshape);

    if nargin < 4 || ~(isscalar(offset) && isnumeric(offset) && offset >= 0)
        offset = 0;
    end

    if ischar(filename)
        filename_ = jrclust.utils.absPath(filename);
        assert(~isempty(filename_), 'could not find file ''%s''', filename);
        fid = fopen(filename_, 'r');
    elseif isnumeric(filename) && isscalar(filename) % fid
        fid = filename;
        assert(ftell(fid) > -1, 'unopened file or invalid fid: %s', fid);
    else
        error('invalid filename or fid: %s', filename);
    end

    fseek(fid, offset, 'bof');

    if nargout > 1
        t = tic;
    end
    % load the data
    if isempty(dshape)
        data = fread(fid, inf, ['*' dtype]);
    elseif ismatrix(dshape) % fread handles this just fine
        data = fread(fid, dshape, ['*' dtype]);
    else
        nSamples = prod(dshape);
        data = fread(fid, nSamples, ['*' dtype]);

        if numel(data) == nSamples
            data = reshape(data, dshape);
        else
            d1 = dshape(1);
            d2 = floor(numel(data)/d1);
            if d2 >= 1
                data = reshape(data, d1, d2);
            else
                data = [];
            end
        end
    end
    fclose(fid);

    if nargout > 1
        loadTime = toc(t);
    end
end
