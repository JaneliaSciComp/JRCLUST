%--------------------------------------------------------------------------
function hFig = figure_new_(vcTag)
    %remove prev tag duplication
    delete_multi_(findobj('Tag', vcTag, 'Type', 'Figure'));

    hFig = figure('Tag', vcTag);
end %func
