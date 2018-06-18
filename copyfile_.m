%--------------------------------------------------------------------------
% 11/5/17 JJJ: added vcDir_from
% 9/26/17 JJJ: multiple targeting copy file. Tested
function copyfile_(csFiles, vcDir_dest, vcDir_from)
    % copyfile_(vcFile, vcDir_dest)
    % copyfile_(csFiles, vcDir_dest)
    % copyfile_(csFiles, csDir_dest)

    if nargin<3, vcDir_from = ''; end
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
        if ~isempty(vcDir_from), vcPath_from_ = fullfile(vcDir_from, vcPath_from_); end
        if exist_dir_(vcPath_from_)
            [vcPath_,~,~] = fileparts(vcPath_from_);
            vcPath_from_ =  sprintf('%s%s*', vcPath_, filesep());
            vcPath_to_ = sprintf('%s%s%s', vcDir_dest, filesep(), dir_filesep_(csFiles{iFile}));
            if ~exist_dir_(vcPath_to_), mkdir(vcPath_to_); end
            %         disp([vcPath_from_, '; ', vcPath_to_]);
        else
            vcPath_to_ = vcDir_dest;
            if ~exist_dir_(vcPath_to_)
                mkdir(vcPath_to_);
                disp(['Created a folder ', vcPath_to_]);
            end
        end
        try
            vcEval1 = sprintf('copyfile ''%s'' ''%s'' f;', vcPath_from_, vcPath_to_);
            eval(vcEval1);
            fprintf('\tCopied ''%s'' to ''%s''\n', vcPath_from_, vcPath_to_);
        catch
            fprintf(2, '\tFailed to copy ''%s''\n', vcPath_from_);
        end
    end
end %func
