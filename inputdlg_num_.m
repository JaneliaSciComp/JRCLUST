%--------------------------------------------------------------------------
function out = inputdlg_num_(vcGuide, vcVal, vrVal)
    csAns = inputdlg_(vcGuide, vcVal, 1, {num2str(vrVal)});
    try
        out = str2double(csAns{1});
    catch
        out = nan;
    end
end %func
