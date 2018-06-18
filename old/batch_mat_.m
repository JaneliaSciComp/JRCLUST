%--------------------------------------------------------------------------
% 10/30/17 JJJ: To be deprecated
function batch_mat_(vcFile_batch_mat, vcCommand)
    % batch process binary file from a template file
    % batch_(myfile_batch.mat, vcCommand)
    %  file must contain: csFiles_bin, csFiles_template
    %    optional: vrDatenum, datenum_start, csFiles_prb

    if ~contains(lower(vcFile_batch_mat), '_batch.mat')
        fprintf(2, 'Must provide _batch.mat file format');
        return;
    end
    if nargin<2, vcCommand = ''; end
    if isempty(vcCommand), vcCommand = 'spikesort'; end

    % Input file format
    S_batch = load(vcFile_batch_mat);
    csFiles_bin = S_batch.csFiles_bin;
    if ischar(csFiles_bin), csFiles_bin = {csFiles_bin}; end
    nFiles = numel(csFiles_bin);
    csFiles_template = get_(S_batch, 'csFiles_template');
    if ischar(csFiles_template), csFiles_template = repmat({csFiles_template}, size(csFiles_bin)); end
    csFiles_prb = get_(S_batch, 'csFiles_prb');
    if isempty(csFiles_prb), csFiles_prb = ''; end
    if ischar(csFiles_prb), csFiles_prb = repmat({csFiles_prb}, size(csFiles_bin)); end
    csFiles_prm = cell(size(csFiles_bin));

    for iFile = 1:nFiles
        try
            vcFile_prm1 = makeprm_template_(csFiles_bin{iFile}, csFiles_template{iFile}, csFiles_prb{iFile});
            fprintf('Created %s\n', vcFile_prm1);
            jrc3(vcCommand, vcFile_prm1);
            csFiles_prm{iFile} = vcFile_prm1;
        catch
            disperr_(sprintf('Failed to process %s', csFiles_bin{iFile}));
        end
    end %for
    S_batch.csFiles_prm = csFiles_prm;
    write_struct_(vcFile_batch_mat, S_batch);
end %func
