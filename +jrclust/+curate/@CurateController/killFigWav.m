function killFigWav(obj, hFig, ~)
%KILLFIGWAV Destroy the main figure, close all other figures and any
%lingering msgboxes.
if ~obj.isEnding
    obj.endSession(); % we'll be back
end

if obj.isEnding
    delete(hFig);
end
end