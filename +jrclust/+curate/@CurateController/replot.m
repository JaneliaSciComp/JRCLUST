function replot(obj)
%REPLOT Refresh FigWav, FigSim, and selected unit
obj.updateFigWav();
obj.updateFigSim();

obj.updateSelect(obj.selected, 1);
end

