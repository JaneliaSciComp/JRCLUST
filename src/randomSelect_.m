%--------------------------------------------------------------------------
function viTime1 = randomSelect_(viTime1, nShow)
    if isempty(viTime1), return; end
    if numel(viTime1) > nShow
        viTime1 = viTime1(randperm(numel(viTime1), nShow));
    end
end %func
