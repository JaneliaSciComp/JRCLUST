%--------------------------------------------------------------------------
function viTime1 = randomSelect_(viTime1, nShow)
    if isempty(viTime1), return; end
    if numel(viTime1) > nShow
        viTime1 = viTime1(ceil(numel(viTime1) .* rand(nShow, 1))); % TW
    end
end %func
