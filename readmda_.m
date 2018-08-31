%--------------------------------------------------------------------------
function A = readmda_(fname)
    % Author: Jeremy Magland, modified by JJJ
    % Jan 2015; Last revision: 15-Feb-2106

    if (strcmp(fname(end-4:end),'.csv')==1)
        A=textread(fname,'','delimiter',',');
        return;
    end

    F=fopen(fname,'r','l');

    % read the first header sample: data type
    try
        code=fread(F,1,'int32');
    catch
        error('Problem reading file: %s',fname);
    end

    % read the second header sample: number of dimensions
    if (code>0)
        num_dims=code;
        code=-1;
    else
        fread(F,1,'int32');
        num_dims=fread(F,1,'int32');
    end;

    % read the length per dimension
    dim_type_str='int32';
    if (num_dims<0)
        num_dims=-num_dims;
        dim_type_str='int64';
    end;

    % read length per dimension
    S = fread(F, num_dims, dim_type_str)';
    N = prod(S);

    switch code
        case -1
        A = fread(F,N*2,'*float');
        A = A(1:2:end) + sqrt(-1) * A(2:2:end);
        case -2, A = fread(F,N,'*uchar');
        case -3, A = fread(F,N,'*float');
        case -4, A = fread(F,N,'*int16');
        case -5, A = fread(F,N,'*int32');
        case -6, A = fread(F,N,'*uint16');
        case -7, A = fread(F,N,'*double');
        case -8, A = fread(F,N,'*uint32');
        otherwise, error('Unsupported data type code: %d',code);
    end
    A = reshape(A, S);
    fclose(F);
end % function
