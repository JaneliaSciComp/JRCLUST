function [hPlot, viY] = multiplot(hPlot, scale, vrX, mrY, viY, fScatter)
    % Create or rescale a multi-line plot
    % hPlot = multiplot([], scale, vrX, mrY) %create a new plot
    % multiplot(hPlot, scale, vrX, mrY) %create a new plot
    % multiplot(hPlot, scale, vrX, mrY, viY) %create a new plot
    % multiplot(hPlot, scale) %rescale. hPlot can be cell or arrays
    % [3D matrix mode]
    % multiplot(hPlot, scale, mrX, trY, miY) %trY: nSamples x nSites x nSpk
    if nargin >= 4
        if isa(mrY, 'gpuArray')
            mrY = gather(mrY);
        end
        if isa(mrY, 'int16')
            mrY = single(mrY);
        end
    end
    if nargin < 6
        fScatter = false;
    end

    if nargin > 2 %create
        shape = size(mrY); % nSamples x nSites x (nSpk)
        if nargin < 5
            viY = 1:shape(2);
        end

        S_plot = struct('scale', scale, 'shape', shape, 'viY', viY, 'fScatter', fScatter);  
        if fScatter % only for points that are not connected
            mrX = vrX(:);
            mrY = mrY(:)/scale + viY(:);        
        elseif ismatrix(mrY)    
            if isempty(vrX)
                vrX = (1:shape(1))';
            end

            if isvector(vrX)
                mrX = repmat(vrX(:), [1, shape(2)]); 
            else
                mrX = vrX;
            end

            if size(mrX,1) > 2
                mrX(end, :) = nan;
            end

            mrY = bsxfun(@plus, mrY/scale, viY(:)');
        else %3D matrix
            if isempty(vrX)
                vrX = bsxfun(@plus, (1:shape(1))', shape(1)*(0:shape(3)-1)); 
            elseif isvector(vrX)
                vrX = repmat(vrX(:), [1, shape(3)]);
            end

            if size(vrX, 1) > 2
                vrX(end,:) = nan;
            end

            if isvector(viY)
                viY = repmat(viY(:), [1, shape(3)]);
            end

            mrX = permute(repmat(vrX, [1, 1, shape(2)]), [1,3,2]); %vrX is mrX        
            mrY = mrY / scale;

            for iSpk = 1:shape(3)
                mrY(:, :, iSpk) = bsxfun(@plus, mrY(:,:,iSpk), viY(:,iSpk)');
            end        
        end

        if isempty(hPlot)
            hPlot = line(mrX(:), mrY(:)); %faster
            set(hPlot, 'UserData', S_plot);
        else
            set(hPlot, 'XData', mrX(:), 'YData', mrY(:), 'UserData', S_plot);
        end
    else %rescale   
        handle_fun_(@rescale_plot_, hPlot, scale);
        viY = [];
    end
end %func

%% LOCAL FUNCTIONS
function rescale_plot_(plotKey, scale)
    % hPlot must have UserData containing scale, shape
    % multi-line plot
    S_plot = get(plotKey, 'UserData');

    % backward copmatible
    if ~isfield(S_plot, 'viY'), S_plot.viY = 1:S_plot.shape(2); end

    % Rescale
    mrY = reshape(get(plotKey, 'YData'), S_plot.shape);
    if isfield(S_plot, 'fScatter')
        fScatter = S_plot.fScatter; 
    else
        fScatter = 0;
    end
    if fScatter
        mrY = (mrY(:) - S_plot.viY(:)) * S_plot.scale; %restore original 
        mrY = mrY(:) / scale + S_plot.viY(:); %convert to plot
    elseif isvector(S_plot.viY)
        mrY = bsxfun(@minus, mrY, S_plot.viY(:)') * S_plot.scale; %restore original 
        mrY = bsxfun(@plus, mrY / scale, S_plot.viY(:)'); %convert to plot
    else
        for iSpk = 1:S_plot.shape(3)
            viY1 = S_plot.viY(:,iSpk)';
            mrY1 = bsxfun(@minus, mrY(:,:,iSpk), viY1) * S_plot.scale;
            mrY(:,:,iSpk) = bsxfun(@plus, mrY1 / scale, viY1); %convert to plot
        end  
    end

    S_plot.scale = scale; %scale is changed
    set(plotKey, 'YData', mrY(:), 'UserData', S_plot);
end

function handle_fun_(hFun, handles, varargin)
    % iterate handle whether it's a cell or matrix
    if isempty(handles), return; end

    if iscell(handles)
        for i=1:numel(handles)
            if numel(handles{i}) == 1
                hFun(handles{i}, varargin{:});
            else
                handle_fun_(hFun, handles{i}, varargin{:});
            end
        end
    else %matrix
        for i=1:numel(handles)
            hFun(handles(i), varargin{:});
        end                
    end
end