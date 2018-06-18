%--------------------------------------------------------------------------
function batch_plot_(vcFile_batch, vcCommand)
    % vcFile_batch: .batch or _batch.mat file format (contains csFiles_prm)
    % Collectively analyze multiple sessions
    % error('not implemented yet');
    if nargin<2, vcCommand=[]; end
    if isempty(vcCommand), vcCommand='skip'; end %spikesort if doesn't exist

    if ~exist(vcFile_batch, 'file'), fprintf(2, 'File does not exist\n'); return; end
    if matchFileExt_(vcFile_batch, '.batch')
        edit(vcFile_batch); %show script
        csFiles_prm = importdata(vcFile_batch);
        % Removing comments that starts with "%"
        func_comment = @(vc)vc(1) == '%';
        viComment = cellfun(@(vc)func_comment(strtrim(vc)), csFiles_prm);
        csFiles_prm(viComment) = [];
    end


    % run the sorting and collect data. quantify the quality
    cS_plot_file = cell(size(csFiles_prm));
    for iFile=1:numel(csFiles_prm)
        try
            vcFile_prm1 = csFiles_prm{iFile};
            jrc3('clear');
            if ~strcmpi(vcCommand, 'skip')
                jrc3(vcCommand, vcFile_prm1);
                S0 = get(0, 'UserData');
            else
                S0 = load_cached_(vcFile_prm1);
            end
            cS_plot_file{iFile} = S_plot_new_(S0);
        catch
            disp(lasterr());
        end
    end %for

    % plot cS_plot_file
    S_plot_show(cS_plot_file); % save to _batch.mat (?)

end %func
