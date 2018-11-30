function mi = rankorderMatrix(meanWf)
    %RANKORDERMATRIX 
    if isrow(meanWf)
        meanWf = meanWf';
    end

    shape = size(meanWf);
    if numel(shape) == 3
        meanWf = meanWf(:); % global order
    end

    mi = zeros(size(meanWf));
    for iCol = 1:size(meanWf, 2)
        col = meanWf(:, iCol);
        pos = find(col > 0);

        if ~isempty(pos)
            [~, argsort] = sort(col(pos), 'ascend');
            mi(pos(argsort), iCol) = 1:numel(pos);
        end

        neg = find(col < 0);
        if ~isempty(neg)
            [~, argsort] = sort(col(neg), 'descend');
            mi(neg(argsort), iCol) = -(1:numel(neg));
        end
    end

    if numel(shape) == 3
        mi = reshape(mi, shape);
    end
end