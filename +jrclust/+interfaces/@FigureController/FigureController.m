classdef (Abstract) FigureController < handle
    %FIGURECONTROLLER A controller driving one or more figures
    properties
    end

    %% LIFECYCLE
    methods
        function obj = FigureController()
            %FIGURECONTROLLER Construct an instance of this class
            %   Detailed explanation goes here
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        menuCheckbox(menuTag, label, fUncheckOthers);
        menuOptions(obj, parent, labels, hFun);
    end
end
