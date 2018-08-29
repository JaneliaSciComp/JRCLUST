%--------------------------------------------------------------------------
% Display list of toolbox and files needed
% 7/26/17 JJJ: Code cleanup and test
function [fList, pList] = disp_dependencies_()
    [fList,pList] = matlab.codetools.requiredFilesAndProducts(mfilename());
    if nargout==0
        disp('Required toolbox:');
        disp({pList.Name}');
        disp('Required files:');
        disp(fList');
    end
end % function
