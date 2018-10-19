%--------------------------------------------------------------------------
% 8/14/17 JJJ: Type Figure added
function hFig = set_fig_(vcTag, S_fig)
    % return figure handle based on the tag
    % hFig = set_fig_(vcTag, S_fig)
    % hFig = set_fig_(hFig, S_fig)
    hFig = [];
    try
        if ischar(vcTag)
            hFig = findobj('Tag', vcTag, 'Type', 'Figure');
        else
            hFig = vcTag;
        end
        set(hFig, 'UserData', S_fig); %figure property
    catch
        disperr_();
    end
end %end
