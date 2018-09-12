%--------------------------------------------------------------------------
% 9/26/17 JJJ: multiple targeting copy file. Tested
function copyfile__(csFiles, vcDir_dest)
    % copyfile_(vcFile, vcDir_dest)
    % copyfile_(csFiles, vcDir_dest)
    % copyfile_(csFiles, csDir_dest)

    % Recursion if cell is used
    if iscell(vcDir_dest)
        csDir_dest = vcDir_dest;
        for iDir = 1:numel(csDir_dest)
            try
                copyfile_(csFiles, csDir_dest{iDir});
            catch
                disperr_();
            end
        end
        return;
    end

    if ischar(csFiles), csFiles = {csFiles}; end
    for iFile=1:numel(csFiles)
        vcPath_from_ = csFiles{iFile};
        if exist(vcPath_from_, 'dir') == 7
            [vcPath_,~,~] = fileparts(vcPath_from_);
            vcPath_from_ =  sprintf('%s%s*', vcPath_, filesep());
            vcPath_to_ = sprintf('%s%s%s%s', vcDir_dest, filesep(), vcPath_, filesep());
            if exist(vcPath_to_, 'dir') ~= 7, mkdir(vcPath_to_); end
            disp([vcPath_from_, '; ', vcPath_to_]);
        else
            vcPath_to_ = vcDir_dest;
            if exist(vcPath_to_, 'dir') ~= 7
                mkdir(vcPath_to_);
                disp(['Created a folder ', vcPath_to_]);
            end
        end
        try
            vcEval1 = sprintf('copyfile ''%s'' ''%s'' f;', vcPath_from_, vcPath_to_);
            eval(vcEval1);
            fprintf('\tCopied %s to %s\n', vcPath_from_, vcPath_to_);
        catch
            fprintf(2, '\tFailed to copy %s\n', vcPath_from_);
        end
    end
end %func
