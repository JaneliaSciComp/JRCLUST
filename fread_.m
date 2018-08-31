%--------------------------------------------------------------------------
function rawTraces = fread_(fidBinary, sampleDims, dataType)
    % Get around fread bug (matlab) where built-in fread resize doesn't work
    try
        if isempty(sampleDims)
            rawTraces = fread(fidBinary, inf, ['*', dataType]);
        else
            if numel(sampleDims) == 1
                sampleDims = [sampleDims, 1];
            end

            rawTraces = fread(fidBinary, prod(sampleDims), ['*', dataType]);

            if numel(rawTraces) == prod(sampleDims)
                rawTraces = reshape(rawTraces, sampleDims);
            else
                dim2 = floor(numel(rawTraces) / sampleDims(1));
                if dim2 >= 1
                    rawTraces = reshape(rawTraces, sampleDims(1), dim2);
                else
                    rawTraces = [];
                end
            end
        end
    catch
        disperr_();
    end
end % function
