%--------------------------------------------------------------------------
function delete_(csFiles)
    if ischar(csFiles), csFiles = {csFiles}; end
    for i=1:numel(csFiles)
        try
            if iscell(csFiles)
                delete(csFiles{i});
            else
                delete(csFiles(i));
            end
        catch
            %         disperr_();
        end
    end
end %func
