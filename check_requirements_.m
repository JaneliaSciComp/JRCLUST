%--------------------------------------------------------------------------
% @TODO: Limit the search es to .m files only
function flag = check_requirements_()
    csToolbox_req = {'distrib_computing_toolbox', 'image_toolbox', 'signal_toolbox', 'statistics_toolbox'};
    [fList, pList] = disp_dependencies_();

    vlToolboxes = cellfun(@(vc)license('test', vc), csToolbox_req);
    vlFiles = cellfun(@(vc)exist(vc, 'file')==2, fList);
    flag = all([vlToolboxes(:); vlFiles(:)]);
end %func
